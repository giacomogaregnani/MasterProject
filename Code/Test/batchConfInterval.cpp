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
    int nMCMC = 50000;

    // DEFINE THE PROBABILISTIC INTEGRATOR
    double sigma = 0.5;

    std::vector<double> hVec = {0.2};
    for (int i = 0; i < 6; i++) {
        hVec.push_back(hVec.back() / 2.0);
    }
    int nMC = 1;
    size_t numberOfReps = 1;

    for (auto h : hVec) {

        size_t expMean = 0;

        std::vector<double> MCMean(numberOfReps);
        std::vector<double> batchMeansEst(numberOfReps);

        std::ofstream thetas;
        std::string thetaFileName = std::string(DATA_PATH) + "/ProbVsDet/" + argv[1] + strTime
                                    + std::to_string(static_cast<int>(h * 1e5)) + "prob.txt";
        thetas.open(thetaFileName, std::ofstream::out | std::ofstream::trunc);

        for (expMean = 0; expMean < numberOfReps; expMean++) {

            // Chains
            std::vector<VectorXd> mcmcPath;
            long int cost = 0;
            std::vector<double> likelihoods = {};

            mcmcPath = MetropolisHastings(odeModel, paramList, sigma,
                                          h, finalTime, data, times,
                                          priorMean, priorVariance, nMC,
                                          nMCMC, varData, &cost, false,
                                          likelihoods, generator);

            MCMean[expMean] = 0.0;
            for (auto it : mcmcPath) {
                MCMean[expMean] += it.dot(it);
            }
            MCMean[expMean] /= nMCMC;

            batchMeansEst[expMean] = computeBatchMeans(mcmcPath, MCMean[expMean]);
        }

        // WRITE RESULTS ON FILE
        for (size_t i = 0; i < numberOfReps; i++) {
            thetas << std::fixed << std::setprecision(30) << batchMeansEst[i] << "\t" << MCMean[i] << "\n";
        }
        thetas.close();

        std::string thetaFileNameDet = std::string(DATA_PATH) + "/ProbVsDet/" + argv[1] + strTime
                                    + std::to_string(static_cast<int>(h * 1e5)) + "det.txt";
        thetas.open(thetaFileNameDet, std::ofstream::out | std::ofstream::trunc);

        for (expMean = 0; expMean < numberOfReps; expMean++) {

            // Chains
            std::vector<VectorXd> mcmcPath;
            long int cost = 0;
            std::vector<double> likelihoods = {};

            mcmcPath = MetropolisHastings(odeModel, paramList, 0.0,
                                          h, finalTime, data, times,
                                          priorMean, priorVariance, nMC,
                                          nMCMC, varData, &cost, false,
                                          likelihoods, generator);

            MCMean[expMean] = 0.0;
            for (auto it : mcmcPath) {
                MCMean[expMean] += it.dot(it);
            }
            MCMean[expMean] /= nMCMC;

            batchMeansEst[expMean] = computeBatchMeans(mcmcPath, MCMean[expMean]);
        }

        // WRITE RESULTS ON FILE
        for (size_t i = 0; i < numberOfReps; i++) {
            thetas << std::fixed << std::setprecision(30) << batchMeansEst[i] << "\t" << MCMean[i] << "\n";
        }
        thetas.close();
    }

    return 0;
}

