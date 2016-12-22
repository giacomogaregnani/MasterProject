#include "mcmcTools.hpp"
#include <fstream>
#include <ctime>

int main()
{
    VectorXd firstGuess(2);
    firstGuess(0) = 1.5;
    firstGuess(1) = -1.0;

    VectorXd fixedFirstGuess = firstGuess;

    std::default_random_engine generator{(unsigned int) time(NULL)};

    double gamma = 0.5;
    int count = 0;

    std::vector<double> rhoBar = {1.1, 1.01, 1.001};

    for (auto it : rhoBar) {
        // Obtain results
        std::vector<std::vector<VectorXd>> MH(6), RAM(6);
        std::vector<VectorXd> MHMix, RAMMix;
        int nMCMC = 100;
        double accRatioMH = 0.0, accRatioRAM = 0.0;
        MatrixXd S = RAMinit(it, 0.4, 2);

        int indexMH[6] = {0, 0, 0, 0, 0, 0};
        int indexRAM[6] = {0, 0, 0, 0, 0, 0};

        // DO IT WITH MH
        bool convergence = false;

        for (size_t j = 0; j < 6; j++) {
            MH[j].push_back(fixedFirstGuess + VectorXd::Random(2));
        }

        while (!convergence) {

            size_t i;
            for (i = 0; i < 6; i++) {
                MH[i] = testMetropolisConv(MH[i].back(), nMCMC, &accRatioMH, gamma,
                                           false, 0.0, S, &(indexMH[i]), generator);
            }

            // MIX
            for (size_t j = 0; j < MH[0].size(); j++) {
                for (size_t u = 0; u < 6; u++) {
                    MHMix.push_back(MH[u][j]);
                }
            }

            // Single variances
            VectorXd tmpMeans;
            VectorXd withinVarianceMean = VectorXd::Zero(2);

            for (size_t j = 0; j < 6; j++) {
                tmpMeans = computeMeans(MH[j]);
                withinVarianceMean += computeVariances(MH[j], tmpMeans);
            }
            withinVarianceMean /= 6;

            // Compute mix variance
            VectorXd mixVariance = VectorXd::Zero(2);
            tmpMeans = computeMeans(MHMix);
            mixVariance = computeVariances(MHMix, tmpMeans);

            // rho-test for convergence
            std::cout << "MH " << indexMH[0] << std::endl;
            convergence = true;
            for (int j = 0; j < 2; j++) {
                double rho = sqrt(mixVariance(j) / withinVarianceMean(j));
                std::cout << rho << " ";
                if (rho > it) {
                    convergence = false;
                }
            }
            std::cout << std::endl;

        }

        // DO IT WITH RAM
        for (size_t j = 0; j < 6; j++) {
            RAM[j].push_back(fixedFirstGuess + VectorXd::Random(2));
        }

        convergence = false;
        while (!convergence) {

            size_t i;
            for (i = 0; i < 6; i++) {
                RAM[i] = testMetropolisConv(RAM[i].back(), nMCMC, &accRatioRAM, gamma,
                                            false, 0.25, S, &(indexRAM[i]), generator);
            }

            // MIX
            for (size_t j = 0; j < RAM[0].size(); j++) {
                for (size_t u = 0; u < 6; u++) {
                    RAMMix.push_back(RAM[u][j]);
                }
            }

            // Single variances
            VectorXd tmpMeans;
            VectorXd withinVarianceMean = VectorXd::Zero(2);

            for (size_t j = 0; j < 6; j++) {
                tmpMeans = computeMeans(RAM[j]);
                withinVarianceMean += computeVariances(RAM[j], tmpMeans);
            }
            withinVarianceMean /= 6;

            // Compute mix variance
            VectorXd mixVariance = VectorXd::Zero(2);
            tmpMeans = computeMeans(RAMMix);
            mixVariance = computeVariances(RAMMix, tmpMeans);

            // rho-test for convergence
            convergence = true;
            std::cout << "RAM " << indexRAM[0] << std::endl;
            for (int j = 0; j < 2; j++) {
                double rho = sqrt(mixVariance(j) / withinVarianceMean(j));
                std::cout << rho << " ";
                if (rho > it) {
                    convergence = false;
                }
            }
            std::cout << std::endl;

        }

        // Write results on file
        std::fstream results, resultsRAM;
        std::string iteration = std::to_string(count++);
        std::string finalpath = std::string(DATA_PATH) + std::string("/resultsMH_") + iteration + std::string(".txt");
        std::string finalpathRAM = std::string(DATA_PATH) + std::string("/resultsRAM_") + iteration + std::string(".txt");

        results.open(finalpath, std::ofstream::out | std::ofstream::trunc);
        resultsRAM.open(finalpathRAM, std::ofstream::out | std::ofstream::trunc);

        for (size_t i = 0; i < MHMix.size(); i++) {
            results << MHMix[i].transpose() << "\n";
        }
        for (size_t i = 0; i < RAMMix.size(); i++) {
            resultsRAM << RAMMix[i].transpose() << "\n";
        }

        results.close();
        resultsRAM.close();
    }

    return 0;
}
