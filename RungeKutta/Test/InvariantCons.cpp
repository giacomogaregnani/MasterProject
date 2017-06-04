#include <RandomTimeStep.hpp>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

double angularMomentum(VectorXd v)
{
    return v(2) * v(1) - v(3) * v(0);
}

double hamKepler(VectorXd v)
{
    return 0.5 * (v(0) * v(0) + v(1) * v(1)) - 1.0 / (sqrt(v(2) * v(2) + v(3) * v(3)));
}

double hamHenon(VectorXd v)
{
    return 0.5 * (v(0) * v(0) + v(1) * v(1)) +
           0.5 * (v(2) * v(2) + v(3) * v(3)) +
           v(2) * v(2) * v(3) - v(3) * v(3) * v(3) / 3.0;
}

#ifndef PI
#define PI 3.14159265359
#endif

int main(int argc, char* argv[])
{
    if (argc < 8) {
        throw std::invalid_argument("Number of inputs must be 6: \n"
                                    "Number of time step reduction\n"
                                    "Max. time step\n"
                                    "Final time\n"
                                    "Number of MC samples\n"
                                    "File name for output\n"
                                    "Write everything (input 'full' "
                                            " do only one time step)"
                                            " or just result \n"
                                    "exponent p\n"
                                    "==========================");
    }


    // ODE
    odeDef ODE;
    ODE.ode = KEPLER;
    setProblem(&ODE);
    std::vector<double> param = ODE.refParam;

    // Numerical method
    Butcher tableau(IMPMID, IMPLICIT, 0);
    double p = std::atof(argv[7]);

    // Integration parameters
    double T = std::atof(argv[3]), h = std::atof(argv[2]);
    unsigned int nMC = static_cast<unsigned int>(std::stoul(argv[4]));
    unsigned int nExp = static_cast<unsigned int>(std::stoul(argv[1]));
    unsigned int N = static_cast<unsigned int>(T / h);

    // Random seed
    std::default_random_engine generator{(unsigned int) time(NULL)};

    std::string outputFilename = std::string(DATA_PATH) + argv[5] + ".txt";
    std::ofstream output;
    output.open(outputFilename, std::ofstream::out | std::ofstream::trunc);

    // Which invariant do I want to look at?
    double (*fInvariant) (VectorXd);
    fInvariant = &hamKepler;

    if (!strcmp(argv[6],"full")) {

        // Initialization
        RungeKuttaRandomH Method(&generator, ODE,
                                 param, tableau, h, p);

        for (unsigned int k = 0; k < nMC; k++) {
            output << std::fixed << std::setprecision(20)
                   << ODE.initialCond.transpose() << "\t"
                   << fInvariant(ODE.initialCond) << std::endl;

            VectorXd solution = ODE.initialCond;

            double numMomentum;

            for (unsigned int i = 0; i < N; i++) {
                solution = Method.oneStep(solution);
                numMomentum = fInvariant(solution);

                if (numMomentum != numMomentum) {
                    output.close();
                    throw std::out_of_range("Solution diverged \n");
                }

                output << solution.transpose() << "\t"
                       << numMomentum << std::endl;
            }
        }

        output.close();

    } else {
        for (unsigned int j = 0; j < nExp; j++) {

            // Initialization
            RungeKuttaRandomH Method(&generator, ODE,
                                      param, tableau, h, p);

            // Output file
            std::cout << std::fixed << std::setprecision(8)
                      << "Computation for h = " << h << std::endl;

            VectorXd solution = ODE.initialCond;
            double refInvariant = fInvariant(ODE.initialCond);
            double errInvariantMax = 0, errInvariantTime = 0;
            double numMomentum;

            for (unsigned int i = 0; i < N; i++) {
                solution = Method.oneStep(solution);
                numMomentum = fInvariant(solution);
                errInvariantTime = std::abs(numMomentum - refInvariant);
                // output << h * i << "\t" << errInvariantTime << std::endl;

                if (errInvariantTime > errInvariantMax) {
                    errInvariantMax = errInvariantTime;
                }

                if (numMomentum != numMomentum) {
                    break;
                }
            }

            output << h << "\t" << errInvariantMax << std::endl;
            std::cout << std::fixed << std::setprecision(20)
                      << "Invariant error = " << errInvariantMax
                      << std::endl;
            h /= 2;
            N *= 2;
        }
        output.close();
    }

    return 0;
}
