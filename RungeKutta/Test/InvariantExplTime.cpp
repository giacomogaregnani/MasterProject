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
    if (argc < 5) {
        throw std::invalid_argument("Number of inputs must be 3: \n"
                                            "Number of time step reduction\n"
                                            "Max. time step\n"
                                            "Number of MC samples\n"
                                            "File name for output\n"
                                            "==========================");
    }

    // ODE
    odeDef ODE;
    ODE.ode = HENHEIL;
    setProblem(&ODE);
    std::vector<double> param = ODE.refParam;

    // Numerical method
    Butcher tableau(IMPMID, IMPLICIT, 0);

    // Integration parameters
    double h = std::atof(argv[2]);
    unsigned int nExp = static_cast<unsigned int>(std::stoul(argv[1]));
    unsigned int nMC = static_cast<unsigned int>(std::stoul(argv[3]));

    // Random seed
    std::default_random_engine generator{(unsigned int) time(NULL)};

    std::string outputFilename = std::string(DATA_PATH) + argv[4] + ".txt";
    std::ofstream output;
    output.open(outputFilename, std::ofstream::out | std::ofstream::trunc);

    // Which invariant do I want to look at?
    double (*fInvariant) (VectorXd);
    fInvariant = &hamHenon;

    for (unsigned int j = 0; j < nExp; j++) {

        std::cout << std::fixed << std::setprecision(6)
                  << "Computation for h = " << h << std::endl
                  << "(" << j << " done, " << nExp - j << " to go)" << std::endl;

        unsigned int k = 0;
        std::vector<double> divTime(nMC);

        unsigned int dead = 0;
        #pragma omp parallel num_threads(30) private(k)
        {
            #pragma omp for
            for (k = 0; k < nMC; k++) {

                // Initialization
                RungeKuttaAddNoise Method(&generator, ODE,
                                          param, tableau, h, 2.0);

                VectorXd solution = ODE.initialCond;
                double numMomentum = 0.0;

                unsigned int i = 0;
                while (numMomentum == numMomentum) {
                    solution = Method.oneStep(solution);
                    numMomentum = fInvariant(solution);
                    i++;
                }
                #pragma omp atomic
                ++dead;

                #pragma omp critical
                std::cout << "Solution diverged. Dead " << dead << std::endl;
                divTime[k] = h * i;
            }
        }

        // Compute mean divergence time
        double meanDivTime = 0.0;
        for (auto it : divTime) {
            meanDivTime += it;
        }
        meanDivTime /= nMC;

        // Output file
        output << h << "\t";

        std::cout << "Solution diverged at time "
                  << meanDivTime << std::endl;
        output << meanDivTime << std::endl;

        h /= 1.2;
    }

    output.close();

    return 0;
}
