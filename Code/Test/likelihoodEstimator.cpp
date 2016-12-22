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

    // TIME STEP AND NUMBER OF MC TRAJECTORIES
    std::vector<double> h = {};
    double maxH = 0.05;
    for (int i = 0; i < 7; i++) {
        h.push_back(maxH);
        maxH /= 2;
    }
    int nMC = 1;

    // NUMBER OF REPETITIONS
    unsigned int nReps = 4;

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
    std::vector<double> paramGuess(nParam);
    for (size_t i = 0; i < nParam; i++) {
        paramGuess[i] = paramList[i] - 0.5;
    }

    // PARAMETERS OF THE CHAIN
    std::vector<int> nMCMCVec = {10};
    for (int i = 0; i < 12; i++) {
        nMCMCVec.push_back(nMCMCVec.back() * 2);
    }
    std::vector<VectorXd> mcmcPath;

    // DEFINE THE PROBABILISTIC INTEGRATOR
    double sigma = 0.5;

    // IS THE PARAMETER POSITIVE?
    bool isPositive = false;
    if (odeModel.ode == BRUSS) {
        isPositive = true;
    }

    // COMPUTE LIKELIHOOD OF THE EXACT SOLUTION
    double exactLikelihood;
    std::vector<VectorXd> exactSol(nData, VectorXd(odeModel.size));
    if (odeModel.ode == POISSON) {
        for (size_t i = 0; i < nData; i++) {
            exactSol[i] = odeModel.exactSol(paramList, times[i]);
        }
    } else {
        std::string exactFileString("refSolFitznagNoNoise");
        std::fstream exactFileIn;
        exactFileIn.open(std::string(DATA_PATH) + "/" + exactFileString + ".txt", std::ios_base::in);

        // Burn the first part of the file
        double burn;
        exactFileIn >> burn;
        exactFileIn >> burn;
        for (size_t i = 0; i < nData; i++) {
            exactFileIn >> burn;
        }

        for (int i = 0; i < nData; i++) {
            for (int j = 0; j < odeModel.size; j++) {
                exactFileIn >> exactSol[i](j);
            }
        }
        exactFileIn.close();
    }
    exactLikelihood = evalLogLikelihood(data, exactSol, odeModel.size, varData);

    // OUTPUT FILEs
    bool printResults = true;
    if (strcmp(argv[1], "debug") == 0) {
        printResults = false;
    }
    std::ofstream results;
    std::string filepath = std::string(DATA_PATH) + "/" + argv[1] + strTime + ".txt";
    std::cout << filepath << std::endl;

    // OUTPUT FILE
    std::ofstream thetaResults;
    std::string thetaFilepath = std::string(DATA_PATH) + "/" + argv[1] + "_thetas_" + strTime + ".txt";
    std::cout << thetaFilepath << std::endl;

    bool allThetas = false;
    if (argc > 2) {
        if (strcmp(argv[2], "thetas") == 0) {
            allThetas = true;
        }
    }

    bool meanThetas = false;
    if (argc > 3) {
        if (strcmp(argv[3], "meantheta") == 0) {
            meanThetas = true;
        }
    }

    double ith = 0.01;
    // COMPUTE THE APPROXIMATED LIKELIHOOD
    for (auto nMCMC = nMCMCVec.end() - 1; nMCMC < nMCMCVec.end(); nMCMC++) {
        std::cout << "Executing timestep : " << *nMCMC << std::endl;
        std::vector<double> allLik(nReps, 0.0);

        std::vector<double> meanTheta(nReps, 0.0);

        size_t i;
        long int cost;

        #pragma omp parallel for num_threads(4) private(i)
        for (i = 0; i < nReps; i++) {
            std::vector<double> likelihoods;

            /* mcmcPath = detMetropolisHastings(odeModel, paramGuess, ith, finalTime,
                                             data, times, priorMean, priorVariance,
                                             nMCMC, varData, &cost, isPositive,
                                             likelihoods, generator);

            if (allThetas) {
                std::ofstream thetas;
                std::string thetaFileName = std::string(DATA_PATH) + "/" + argv[1] + strTime
                                            + std::to_string(static_cast<int>(ith * 1e6)) + "_det.txt";
                thetas.open(thetaFileName, std::ofstream::out | std::ofstream::trunc);
                for (auto it : mcmcPath) {
                    thetas << it.transpose() << "\n";
                }
                thetas.close();
            } */

            printf("%d\n", i);

            mcmcPath = MetropolisHastings(odeModel, paramGuess, sigma, ith, finalTime,
                                          data, times, priorMean, priorVariance, nMC,
                                          *nMCMC, varData, &cost, isPositive,
                                          likelihoods, generator);
            for (auto it : likelihoods) {
                allLik[i] += it;
            }

            allLik[i] /= *nMCMC;

            if (allThetas) {
                std::ofstream thetas;
                std::string thetaFileName = std::string(DATA_PATH) + "/" + argv[1] + strTime
                                            + std::to_string(*nMCMC) + ".txt";
                thetas.open(thetaFileName, std::ofstream::out | std::ofstream::trunc);
                for (auto it : mcmcPath) {
                    thetas << it.dot(it) << "\n";
                }
                thetas.close();
            }

            if (meanThetas) {
                for (auto it : mcmcPath) {
                    meanTheta[i] += it.dot(it);
                }
                meanTheta[i] /= *nMCMC;
            }

        }

        // Compute variance and bias of the likelihood estimator
        double meanLikelihood = 0.0, varianceLikelihood = 0.0, tmp;
        for (auto it : allLik) {
            meanLikelihood += it;
        }
        meanLikelihood /= nReps;
        double biasSqd = (meanLikelihood - exactLikelihood);
        biasSqd *= biasSqd;

        for (auto it : allLik) {
            tmp = it - meanLikelihood;
            varianceLikelihood += tmp * tmp;
        }
        varianceLikelihood /= nReps;

        if (printResults) {
            results.open(filepath, std::ofstream::out | std::ofstream::app);
            results << std::fixed << std::setprecision(30) << *nMCMC << " " << varianceLikelihood << " " << biasSqd
                    << "\n";
            //results << std::fixed << std::setprecision(30) << ith << " " << varianceLikelihood << " " << biasSqd
            //        << "\n";
            results.close();
        }

        // Compute variance and mean of the MCMC estimator
        double meanMCMC = 0.0, varMCMC = 0.0;
        for (auto it : meanTheta) {
            meanMCMC += it;
        }
        meanMCMC /= nReps;

        if (printResults) thetaResults.open(thetaFilepath, std::ofstream::out | std::ios_base::app);

        for (auto it : meanTheta) {
            tmp = it - meanMCMC;
            varMCMC += tmp * tmp;
        }
        varMCMC /= nReps;

        thetaResults << std::fixed << std::setprecision(30) << *nMCMC << " " << meanMCMC << " " << varMCMC << "\n";
        // thetaResults << std::fixed << std::setprecision(30) << ith << " " << meanMCMC << " " << varMCMC << "\n";
        if (printResults) thetaResults.close();
    }

    if (printResults) {
        results.close();
    }

    return 0;
}