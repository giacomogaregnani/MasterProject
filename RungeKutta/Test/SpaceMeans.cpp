#include <RandomTimeStep.hpp>
#include <UtilitiesGG.hpp>

#include <iomanip>

double phi(VectorXd x)
{
    return x.dot(x);
}

int main(int argc, char* argv[]) {
    // Set problem and solver features
    odeDef ODE;
    ODE.ode = LORENZ;
    setProblem(&ODE);
    std::vector<double> param = ODE.refParam;
    Butcher tableau(EXPTRAPEZ, EXPLICIT, 0);
    double T = 100.0, h = 0.05, p = 1.5;
    unsigned int nMC = 10000, nExp =  3, N;

    // Initialize the random generator
    std::default_random_engine generator{(unsigned int) time(NULL)};

    for (size_t j = 0; j < nExp; j++) {

        std::cout << "========= INFO ========="
                  << std::fixed << std::setprecision(4)
                  << std::endl << "timestep = " << h
                  << std::endl
                  << "========================"
                  << std::endl;

        N = static_cast<unsigned int>(T / h);

        // Information file (useful for plotting in MATLAB)
        std::string fileNameInfo = DATA_PATH + std::string("ResultsSpaceVsTime/infoFile")
                                   + std::to_string(j + 1) + std::string(".txt");
        std::ofstream fileInfo(fileNameInfo, std::ofstream::out | std::ofstream::trunc);
        fileInfo << nMC << "\n" << T << "\n" << h << "\n";
        fileInfo.close();

        std::cout << "========= SPACE AVERAGE COMPUTATION (prob) =========" << std::endl;

        // Initialize nMC solvers, solutions, and space average
        std::vector<RungeKuttaRandomH> solversVector(nMC, RungeKuttaRandomH(&generator, ODE,
                                                                            param, tableau, h, p));
        std::vector<VectorXd> solutionsVector(nMC, VectorXd(ODE.size));
        double spaceAverage;
        for (size_t k = 0; k < nMC; k++) {
            solutionsVector[k] = ODE.initialCond;
        }

        // Open output file
        std::string fileNameSpace = DATA_PATH + std::string("ResultsSpaceVsTime/SolutionSpace")
                                    + std::to_string(static_cast<int>(h * 10000)) + std::string(".txt");
        std::ofstream solStepFile(fileNameSpace, std::ofstream::out | std::ofstream::trunc);

        // An index
        unsigned int k;

        // Loop over time
        for (size_t i = 0; i < N; i++) {

            spaceAverage = 0.0;
            // Do one step for all the solvers and update space averages
            #pragma omp parallel for num_threads(30) private(k) reduction(+:spaceAverage)
            for (k = 0; k < nMC; k++) {
                solutionsVector[k] = solversVector[k].oneStep(solutionsVector[k]);
                spaceAverage += phi(solutionsVector[k]);
            }
            spaceAverage /= nMC;

            solStepFile << spaceAverage << "\n";

            if (i % 100 == 0) {
                std::cout << std::fixed << std::setprecision(1)
                          << "Time = " << i * h
                          << "s, out of " << T << "s" << std::endl;
            }

        }

        // Close file
        solStepFile.close();

        h /= 2;
    }


    return 0;
}