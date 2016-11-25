#include <cmath>
#include <vector>
#include <random>
#include <Solver.hpp>
#include <Eigen/Dense>
#include <ctime>
#include <iostream>
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
    double maxH = 0.1;
    for (int i = 0; i < 14; i++) {
        h.push_back(maxH);
        maxH /= 2;
    }
    int nMC = 10;

    // NUMBER OF REPETITIONS
    unsigned int nReps = 40;

    // PROBLEM DATA
    odeDef odeModel;
    odeModel.ode = POISSON;
    setProblem(&odeModel);
    std::vector<double> paramList = odeModel.refParam;

    // DATA ACQUISITION FROM REFSOL
    std::string refFile("refSolPoisson");
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
        paramGuess[i] = paramList[i];
    }

    // PARAMETERS OF THE CHAIN
    std::vector<VectorXd> mcmcPath;
    int nMCMC = 1;

    // DEFINE THE PROBABILISTIC INTEGRATOR
    double sigma = 0.5;

    // IS THE PARAMETER POSITIVE?
    bool isPositive = false;
    if (odeModel.ode == BRUSS) {
        isPositive = true;
    }

    // COST (UNUSED)
    long int cost;

    // IN THE POISSON CASE THE EXACT SOLUTION IS AVAILABLE, LET'S COMPUTE THE EXACT LIKELIHOOD
    if (odeModel.ode == POISSON) {
        double exactLikelihood;
        std::vector<VectorXd> exactSol(nData, VectorXd(odeModel.size));
        for (size_t i = 0; i < nData; i++) {
            exactSol[i] = odeModel.exactSol(paramList, times[i]);
        }
        exactLikelihood = evalLogLikelihood(data, exactSol, odeModel.size, varData);
        std::ofstream exactFile;
        std::string exactpath = std::string(DATA_PATH) + "/" + argv[1] + strTime
                                + "exact.txt";
        exactFile.open(exactpath, std::ofstream::out | std::ofstream::trunc);
        exactFile << std::fixed << std::setprecision(16) << exactLikelihood;
        exactFile.close();
    }

    for (auto ith : h) {

        std::cout << "Executing timestep : " << ith << std::endl;

        // OUTPUT FILE
        bool printResults = true;
        if (strcmp(argv[1], "debug") == 0) {
            printResults = false;
        }
        std::ofstream results;
        std::string filepath = std::string(DATA_PATH) + "/" + argv[1] + strTime
                               + "_h_"
                               + std::to_string(static_cast<int>(ith*1e9)) + ".txt";
        if (printResults) results.open(filepath, std::ofstream::out | std::ofstream::trunc);
        std::cout << filepath << std::endl;

        // COMPUTE THE APPROXIMATED LIKELIHOOD
        std::vector<double> allLik(nReps);
        unsigned int i;
        #pragma omp parallel for num_threads(20) private(i)
        for (i = 0; i < nReps; i++) {
            std::vector<double> likelihoods = {};
            mcmcPath = MetropolisHastings(odeModel, paramGuess, sigma, ith, finalTime,
                                          data, times, priorMean, priorVariance, nMC,
                                          nMCMC, varData, &cost, isPositive,
                                          likelihoods, generator);
            allLik[i] = likelihoods[0];
        }

        for (auto it : allLik) {
            results << std::fixed << std::setprecision(16) << it << "\n";
        }

        if (printResults) results.close();

    }

    return 0;
}