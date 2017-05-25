#include <RandomTimeStep.hpp>
#include <UtilitiesGG.hpp>

#include <iomanip>

int main ()
{
    // Set problem and solver features
    odeDef ODE;
    ODE.ode = LORENZ;
    setProblem(&ODE);
    Butcher tableau(RK4, EXPLICIT, 0);
    double T = 100.0, h = 0.01;
    unsigned int nMC = 10000, N = static_cast<unsigned int>(T / h);

    // Random generator
    std::default_random_engine generator{(unsigned int) time(NULL)};
    double stdDev = 0.1;

    // Information file (useful for plotting in MATLAB)
    std::string fileNameInfo = DATA_PATH + std::string("ResultsICLorenz/infoFile.txt");
    std::ofstream fileInfo(fileNameInfo, std::ofstream::out | std::ofstream::trunc);
    fileInfo << T << "\n" << h << "\n";
    fileInfo.close();

    // Report progress of for loop
    auto stepSize   = 1ul;
    size_t stepsCompleted = 0;

    // Save final values in a vector
    std::vector<VectorXd> finalValues(nMC, VectorXd(ODE.size));

    unsigned int k;
    // Monte Carlo loop
    #pragma omp parallel
    {
        size_t localCount = 0;

        #pragma omp for private(k)
        for (k = 0; k < nMC; k++) {

            RungeKutta Method(ODE, ODE.refParam, tableau);
            VectorXd solution = ODE.initialCond +
                                normalZeroMeanRandVec(ODE.size, &generator, stdDev);

            // time integration
            for (unsigned int i = 0; i < N; i++) {
                solution = Method.oneStep(solution, h);
            }

            finalValues[k] = solution;

            // Report progress of for loop
            if (localCount++ % stepSize == stepSize - 1) {
                #pragma omp atomic
                ++stepsCompleted;

                #pragma omp critical
                std::cout << "Progress: " << stepsCompleted << " of " << nMC
                          << " (" << std::fixed << std::setprecision(1)
                          << (100.0 * stepsCompleted / nMC) << "%)\n";
            }
        }
    }

    // Output file
    std::string outputFileName = DATA_PATH + std::string("ResultsICLorenz/ICLorenz.txt");
    std::ofstream outputFile(outputFileName, std::ofstream::out | std::ofstream::trunc);
    for (auto it : finalValues) {
        outputFile << it.transpose() << "\n";
    }
    outputFile.close();

    return 0;
}