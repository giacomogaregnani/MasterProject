#include <fstream>
#include <iostream>
#include <iomanip>
#include "generateObservations.hpp"
#include "computeHomogeneous.hpp"
#include <MCMC.hpp>
#include "../matplotlib-cpp-master/matplotlibcpp.h"

namespace plt = matplotlibcpp;

// Remark: epsilon = p(0)

double gradV0(double x)
{
    return x;
}

double V1(double x)
{
    return std::cos(x);
}

double gradV1(double x)
{
    return -std::sin(x);
}

double multiDrift(double x, VectorXd& p)
{
    return -1.0 * (std::exp(p(1)) * gradV0(x) + 1 / p(0) * gradV1(x / p(0)));
}

double homoDrift(double x, VectorXd& p)
{
    return -1.0 * std::exp(p(1)) * gradV0(x);
}

double diffusion(double x, VectorXd &p)
{
    return std::sqrt(2.0 * std::exp(p(2)));
}

int main(int argc, char* argv[])
{
    oneDimSde sde;
    sde.drift = &multiDrift;
    sde.diffusion = &diffusion;
    oneDimSde sdeHomo;
    sdeHomo.drift = &homoDrift;
    sdeHomo.diffusion = &diffusion;

    double T = 10;
    unsigned int N = 1000;

    VectorXd tmpParam(3);
    tmpParam(0) = 0.1;  // Epsilon
    tmpParam(1) = 1.0;   // True multiscale alpha
    tmpParam(2) = 0.5;   // True multiscale betainv

    // std::ofstream output(DATA_PATH + std::string("MultiMulti.txt"), std::ofstream::out | std::ofstream::trunc);
    // std::ofstream outputSol(DATA_PATH + std::string("MultiMultiSol.txt"), std::ofstream::out | std::ofstream::trunc);

    std::ofstream output(DATA_PATH + std::string("MultiHomo.txt"), std::ofstream::out | std::ofstream::trunc);
    std::ofstream outputSol(DATA_PATH + std::string("MultiHomoSol.txt"), std::ofstream::out | std::ofstream::trunc);

    // Compute coefficients of the homogenised equation
    std::vector<double> homCoeffs = computeHomCoeffs(tmpParam, (2.0 * M_PI), &V1);
    output << homCoeffs[0] << "\t" << homCoeffs[1] << std::endl;

    // ============= Transform for positiveness ============= //
    VectorXd param = tmpParam;
    param(1) = std::log(param(1));
    param(2) = std::log(param(2));
    // ====================================================== //

    // Initialize structures for the inverse problem
    unsigned long M = 40, nMCMC = 1001;
    double noise = 1e-4;
    double IC = 5.0;
    std::random_device dev;
    std::default_random_engine noiseSeed{dev()};
    std::default_random_engine proposalSeed{dev()};
    std::default_random_engine acceptanceSeed{dev()};
    std::normal_distribution<double> noiseDistribution(0.0, noise);
    bool IS = true;

    // Generate and perturb observations
    auto x = generateObservations1D(sde, IC, param, T, N, 0);
    for (auto const &itSol : x)
        outputSol << std::fixed << std::setprecision(5) << itSol << "\t";
    outputSol << std::endl;
    outputSol << x[0] << "\t";
    for (unsigned long i = 1; i < x.size(); i++) {
        x[i] += noiseDistribution(noiseSeed);
        outputSol << x[i] << "\t";
    }
    outputSol << std::endl;

    // Compute a homogeneous path with the same Brownian for comparison
    VectorXd homoParam(3);
    homoParam(0) = 0.1;
    homoParam(1) = std::log(homCoeffs[0]);
    homoParam(2) = std::log(homCoeffs[1]);
    auto xHom = generateObservations1D(sdeHomo, IC, homoParam, T, N, 0);
    for (auto const &itSol : xHom) {
        outputSol << std::fixed << std::setprecision(10) << itSol << "\t";
    }
    outputSol << std::endl;

    // Compute the modeling error statistics and rescale the observations
    bool doModErr = false;
    std::vector<double> means, stdDevs;
    std::vector<double> rescaledObs = x;
    std::vector<double> timeVec(N+1);
    for (unsigned int i = 0; i < N+1; i++) {
        timeVec[i] = T/N * i;
    }
    if (doModErr) {
        std::cout << "Computing modeling error..." << std::endl;
        VectorXd priorMean(param.size());
        priorMean << param(0), param(1), param(2);
        VectorXd priorStdDev(param.size());
        priorStdDev << 0.0, 1.0, 1.0;
        ModErr modErr(sdeHomo, sde, &V1, IC, priorMean, priorStdDev, T, N, x, noise);
        unsigned int nMC = 150;
        unsigned int nParam = 40;
        modErr.computePF(nParam, nMC);
        modErr.getStats(means, stdDevs);
        std::cout << "Computed modeling error" << std::endl;
        for (unsigned int i = 0; i < N+1; i++) {
            rescaledObs[i] = x[i] - means[i];
        }
        plt::named_plot("xe", timeVec, x, "b");
        plt::named_plot("x0", timeVec, xHom, "r");
        plt::named_plot("xt", timeVec, rescaledObs, "k");
        plt::legend();
        plt::show();

        for (auto const &itSol : rescaledObs) {
            outputSol << std::fixed << std::setprecision(10) << itSol << "\t";
        }
        outputSol << std::endl;
        for (auto const &it : stdDevs) {
            outputSol << std::fixed << std::setprecision(10) << it << "\t";
        }
        outputSol << std::endl;
    }

    // Initial parameter guess
    VectorXd initGuess = VectorXd::Zero(param.size());
    initGuess(0) = param(0);

    // Inverse problem
    std::shared_ptr<Posterior> posterior;
    posterior = std::make_shared<PFPosterior>(rescaledObs, T, IC, 1, noise, sde, param(0), M, IS); //, stdDevs);
    // posterior = std::make_shared<PFPosteriorHom>(rescaledObs, T, IC, 1, noise, sdeHomo, &V1, param(0), M, IS); //, stdDevs);
    std::shared_ptr<Proposals> proposal;
    std::vector<double> factors = {1.0, 2.0, 1.0};
    proposal = std::make_shared<Proposals>(5e-2, factors);
    MCMC mcmc(initGuess, proposal, posterior, nMCMC);
    auto sample = mcmc.compute(&proposalSeed, &acceptanceSeed);
    for (auto const &itSample : sample)
        output << itSample.transpose() << std::endl;

    output.close();
    outputSol.close();
    return 0;
}