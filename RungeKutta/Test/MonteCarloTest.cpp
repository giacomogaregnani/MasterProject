#include <RandomTimeStep.hpp>
#include <iostream>
#include <fstream>
#include <iomanip>

double phi(VectorXd v)
{
    return v.dot(v);
}

int main(int argc, char* argv[])
{
    if (argc < 4) {
        throw std::invalid_argument("Number of inputs must be 3: \n"
                                            "Output file name\n"
                                            "Reference solution file name\n"
                                            "p\n"
                                            "==========================");
    }

    odeDef ODE;
    ODE.ode = FITZNAG;
    setProblem(&ODE);

    std::vector<double> param = {0.2, 0.2, 3.0};

    Butcher tableau(EXPTRAPEZ, EXPLICIT, 0);

    double T = 10.0, h = 0.1;

    std::default_random_engine generator{(unsigned int) time(NULL)};

    // Output file
    std::string outputFilename = std::string(DATA_PATH) + argv[1] + ".txt";
    std::ofstream output;
    output.open(outputFilename, std::ofstream::out | std::ofstream::trunc);

    // Reference solution
    VectorXd refSolution(ODE.size);
    if (ODE.ode == TEST1D) {
        refSolution(0) = exp(param[0] * T);
    } else {
        std::string refFilename = std::string(DATA_PATH) + argv[2] + ".txt";
        std::fstream reference;
        reference.open(refFilename, std::ofstream::in);
        for (int i = 0; i < ODE.size; i++) {
            reference >> refSolution(i);
        }
    }
    double refValue = phi(refSolution);

    // Monte Carlo parameters
    unsigned int nTimeSteps = 6, nRepetitions = 300, nMC = 1, k;
    std::vector<double> MSEReps(nRepetitions), MSE(nTimeSteps);

    // Loop on the time steps
    for (size_t i = 0; i < nTimeSteps; i++) {

        size_t N = static_cast<size_t>(T / h);
        std::cout << "time step = " << h << std::endl;

        // Repeat the computation of the estimator
        #pragma omp parallel for num_threads(30) private(k)
        for (k = 0; k < nRepetitions; k++) {

            RungeKuttaRandomH Method(&generator, ODE, param, tableau, h, std::atoi(argv[3]));
            double MCEstimator = 0.0;

            // Monte Carlo loop
            for (size_t j = 0; j < nMC; j++) {

                VectorXd solution = ODE.initialCond;

                // Time loop
                for (size_t t = 0; t < N; t++) {
                    solution = Method.oneStep(solution);
                }

                MCEstimator += phi(solution);
            }

            MCEstimator /= nMC;
            MSEReps[k] = (MCEstimator - refValue) * (MCEstimator - refValue);
        }

        MSE[i] = 0.0;
        for (auto it : MSEReps) {
            MSE[i] += it;
        }
        MSE[i] /= nRepetitions;

        output << std::fixed << std::setprecision(30)
               << h << "\t" << MSE[i] << "\n";

        h /= 2;
    }

    double meanOrder = 0.0;
    for (size_t j = 0; j < nTimeSteps; j++) {
        if (j > 0) {
            std::cout << "MSE = " << MSE[j] << "\t";
            std::cout << "order = " << std::log(MSE[j - 1] / MSE[j]) / std::log(2.0) << std::endl;
            meanOrder += std::log(MSE[j - 1] / MSE[j]) / std::log(2.0);
        } else {
            std::cout << "MSE = " << MSE[j] << std::endl;
        }
    }
    meanOrder /= nTimeSteps - 1;
    std::cout << "Mean order of convergence = " << meanOrder << std::endl;


    output.close();

    return 0;
}