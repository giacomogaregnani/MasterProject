#include <RandomTimeStep.hpp>
#include <UtilitiesGG.hpp>

#include <iomanip>

int main(int argc, char* argv[])
{
    // Write to file all values or only the last
    bool fullWrite = false;
    if (argc == 1) {
        fullWrite = true;
    }

    // Set problem and solver features
    odeDef ODE;
    ODE.ode = PEROX;
    setProblem(&ODE);
    std::vector<double> param = ODE.refParam;
    Butcher tableau(RKC, STABEXP, 150);
    double T = 200.0, h = 0.05, p = 1.0;
    unsigned int nMC = 30, N = static_cast<unsigned int>(T / h);

    // Initialize the random generator
    std::default_random_engine generator{(unsigned int) time(NULL)};

    // Delete all files in temporary results folder
    system("exec rm -r " DATA_PATH "/ResultsAdd/*");
    system("exec rm -r " DATA_PATH "/ResultsStep/*");

    // Information file (useful for plotting in MATLAB)
    std::string fileNameInfo = DATA_PATH + std::string("ResultsAdd/infoFile.txt");
    std::ofstream fileInfo(fileNameInfo, std::ofstream::out | std::ofstream::trunc);
    fileInfo << nMC << "\n" << T << "\n" << h << "\n" << fullWrite;
    fileInfo.close();

    // Try with a lot of trajectories
    unsigned int k;

    // Report progress of for loop
    auto stepSize   = 1ul;
    size_t stepsCompleted = 0;

    // The if-else statement has only a difference of file writing
    if (fullWrite) {
        #pragma omp parallel
        {
            size_t localCount = 0;

            #pragma omp for private(k)
            for (k = 0; k < nMC; k++) {

                std::string fileNameAdd = DATA_PATH + std::string("ResultsAdd/SolutionAdd") +
                                          std::to_string(k) + std::string(".txt");
                std::string fileNameStep = DATA_PATH + std::string("ResultsStep/SolutionStep") +
                                           std::to_string(k) + std::string(".txt");
                std::ofstream solAddFile(fileNameAdd, std::ofstream::out | std::ofstream::trunc);
                std::ofstream solStepFile(fileNameStep, std::ofstream::out | std::ofstream::trunc);


                // Additive noise method
                RungeKuttaAddNoise MethodAdd(&generator, ODE, param, tableau, h, p);
                VectorXd solutionAdd = ODE.initialCond;

                // Random time-stepping method
                RungeKuttaRandomH MethodStep(&generator, ODE, param, tableau, h, p + 0.5);
                VectorXd solutionStep = ODE.initialCond;

                // time integration
                for (unsigned int i = 0; i < N; i++) {
                    solAddFile << solutionAdd.transpose() << "\n";
                    solStepFile << solutionStep.transpose() << "\n";

                    solutionAdd = MethodAdd.oneStep(solutionAdd);
                    solutionStep = MethodStep.oneStep(solutionStep);
                }

                solAddFile << solutionAdd.transpose() << "\n";
                solStepFile << solutionStep.transpose() << "\n";

                // Close files
                solAddFile.close();
                solStepFile.close();

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
    } else {

        std::string fileNameAdd = DATA_PATH + std::string("ResultsAdd/SolutionAdd.txt");
        std::string fileNameStep = DATA_PATH + std::string("ResultsStep/SolutionStep.txt");
        std::ofstream solAddFile(fileNameAdd, std::ofstream::out | std::ofstream::trunc);
        std::ofstream solStepFile(fileNameStep, std::ofstream::out | std::ofstream::trunc);

        std::vector<VectorXd> finalTimeResultsAdd(nMC, VectorXd(ODE.size));
        std::vector<VectorXd> finalTimeResultsStep(nMC, VectorXd(ODE.size));


        #pragma omp parallel
        {
            size_t localCount = 0;

            #pragma omp for private(k)
            for (k = 0; k < nMC; k++) {

                // Additive noise method
                RungeKuttaAddNoise MethodAdd(&generator, ODE, param, tableau, h, p);
                VectorXd solutionAdd = ODE.initialCond;

                // Random time-stepping method
                RungeKuttaRandomH MethodStep(&generator, ODE, param, tableau, h, p + 0.5);
                VectorXd solutionStep = ODE.initialCond;

                // time integration
                for (unsigned int i = 0; i < N; i++) {
                    solutionAdd = MethodAdd.oneStep(solutionAdd);
                    solutionStep = MethodStep.oneStep(solutionStep);
                }

                finalTimeResultsAdd[k] = solutionAdd;
                finalTimeResultsStep[k] = solutionStep;

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

        for (size_t i = 0; i < nMC; i++) {
            solAddFile << finalTimeResultsAdd[i].transpose() << "\n";
            solStepFile << finalTimeResultsStep[i].transpose() << "\n";
        }

        // Close files
        solAddFile.close();
        solStepFile.close();
    }


    return 0;
}