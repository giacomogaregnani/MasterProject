#include <cmath>
#include <vector>
#include <random>
#include <Solver.hpp>
#include <fstream>
#include "problems.hpp"
#include "mcmcTools.hpp"
#include <chrono>
#include <iomanip>

using namespace Eigen;

#define PI 3.1415926535897

int main(int argc, char* argv[])
{
    // NUMBER OF PARALLEL MCMC CHAINS
    unsigned int nReps = 6;

    // PROBLEM DATA
    odeDef odeModel;
    odeModel.ode = FITZNAG;
    setProblem(&odeModel);
    std::vector<double> paramList = odeModel.refParam;

    // DATA ACQUISITION FROM REFSOL
    std::string refFile("refSolFitznag");
    std::fstream refSolution;
    refSolution.open(std::string(DATA_PATH) + "/" + refFile + ".txt", std::ios_base::in);

    double finalTime;
    int nData;
    refSolution >> finalTime;
    refSolution >> nData;

    std::vector<double> times(static_cast<unsigned int>(nData));
    std::vector<VectorXd> data(static_cast<unsigned int>(nData), VectorXd(odeModel.size));
    for (int i = 0; i < nData; i++) {
        refSolution >> times[i];
    }
    for (int i = 0; i < nData; i++) {
        for (int j = 0; j < odeModel.size; j++) {
            refSolution >> data[i](j);
        }
    }
    refSolution.close();

    // DATA UNCERTAINTY (COHERENT WITH REFSOL)
    double varData = 1e-2;
    std::default_random_engine generator{(unsigned int) time(NULL)};

    // GET THE DATE TO PUT IT IN THE NAME OF OUTPUT FILE
    time_t now;
    now = time(NULL);
    char charDate[18];
    if (now != -1) {
        strftime(charDate, 18, "%d_%m_%Y_%I_%M_", gmtime(&now));
    }
    std::string strTime(charDate);

    // PRIOR MEAN AND VARIANCE
    size_t nParam = paramList.size();
    VectorXd priorMean(nParam), priorVariance(nParam);
    for (size_t i = 0; i < nParam; i++) {
        priorMean(i) = paramList[i];
        priorVariance(i) = 1.0;
    }

    // INITIAL MCMC GUESS
    VectorXd paramGuess(nParam);
    for (size_t i = 0; i < nParam; i++) {
        paramGuess(i) = paramList[i];
    }

    // PARAMETERS OF THE CHAIN
    std::vector<int> nMCMCVec = {10};
    for (int i = 0; i < 12; i++) {
        nMCMCVec.push_back(nMCMCVec.back() * 2);
    }
    int nMCMC = 500;

    // DEFINE THE PROBABILISTIC INTEGRATOR
    double sigma = 0.5;

    // IS THE PARAMETER POSITIVE?
    bool isPositive = false;
    if (odeModel.ode == BRUSS) {
        isPositive = true;
    }

    // Define the implicit method
    Butcher butcher(GAUSS, 2);

    double h = 0.1;
    int nMC = 1;

    int numberOfReps = 10;

    for (int k = 0; k < 8; k++) {

        h /= 2;
        // numberOfReps *= 2;

        for (int expMean = 0; expMean < numberOfReps; expMean++) {

            std::cout << "h = " << h << std::endl;

            // Chains
            std::vector<std::vector<VectorXd>> mcmcPaths(nReps);
            for (size_t j = 0; j < nReps; j++) {
                mcmcPaths[j].push_back(paramGuess + 0.1 * VectorXd::Random(nParam));
            }
            std::vector<VectorXd> mcmcMix;

            // RAM
            double gamma = 0.01;
            double desiredAlpha = 0.25;
            std::vector<MatrixXd> S(nReps);
            for (size_t j = 0; j < nReps; j++) {
                S[j] = RAMinit(gamma, desiredAlpha, nParam);
            }

            bool convergence = false;

            int count = 0;
            while (!convergence) {
                size_t i;

                #pragma omp parallel for num_threads(6) private(i)
                for (i = 0; i < nReps; i++) {
                    std::vector<double> stdInGuess(nParam);
                    for (size_t j = 0; j < nParam; j++) {
                        stdInGuess[j] = (mcmcPaths[i].back())(j);
                    }


                    std::vector<double> likelihoods;
                    printf("%zu\n", i);
                    StoppedMetropolis(mcmcPaths[i], odeModel, stdInGuess, sigma, h, finalTime,
                                      data, times, priorMean, priorVariance, nMC,
                                      nMCMC, varData, isPositive,
                                      likelihoods, generator, S[i], count * nMCMC);
                    /* impStoppedMetropolis(mcmcPaths[i], odeModel, stdInGuess, sigma, h, finalTime,
                                         data, times, priorMean, priorVariance, nMC,
                                         nMCMC, varData, isPositive,
                                         likelihoods, generator, S[i], count * nMCMC, butcher); */
                }
                count++;
                #pragma omp barrier

                // Mix (interleaving)
                mcmcMix.clear();

                for (size_t j = 0; j < mcmcPaths[0].size(); j++) {
                    for (size_t u = 0; u < nReps; u++) {
                        mcmcMix.push_back(mcmcPaths[u][j]);
                    }
                }

                // Compute single variances
                VectorXd tmpMeans;
                VectorXd withinVarianceMean = VectorXd::Zero(nParam);

                for (size_t j = 0; j < nReps; j++) {
                    tmpMeans = computeMeans(mcmcPaths[j]);
                    withinVarianceMean += computeVariances(mcmcPaths[j], tmpMeans);
                }
                withinVarianceMean /= nReps;

                // Compute mix variance
                VectorXd mixVariance = VectorXd::Zero(nParam);
                tmpMeans = computeMeans(mcmcMix);
                mixVariance = computeVariances(mcmcMix, tmpMeans);

                // rho-test for convergence
                convergence = true;
                for (int j = 0; j < nParam; j++) {
                    double rho = sqrt(mixVariance(j) / withinVarianceMean(j));
                    std::cout << rho << " ";
                    if (rho > 1.03) {
                        convergence = false;
                    }
                }
                std::cout << std::endl;
            }

            // WRITE RESULTS ON FILE
            std::ofstream thetas;
            std::string thetaFileName = std::string(DATA_PATH) + "/MeanExp/" + argv[1] + strTime
                                        + std::to_string(static_cast<int>(h * 1e5)) + "_"
                                        + std::to_string(expMean) + ".txt";
            thetas.open(thetaFileName, std::ofstream::out | std::ofstream::trunc);
            for (auto it : mcmcMix) {
                thetas << it.transpose() << "\n";
            }
            thetas.close();
        }
    }

    return 0;
}
