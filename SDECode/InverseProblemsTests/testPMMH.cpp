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
    oneDimSde sde{&multiDrift, &diffusion};
    oneDimSde sdeHomo{&homoDrift, &diffusion};

    double T = 200;
    unsigned int N = 40000;

    VectorXd tmpParam(3);
    tmpParam(0) = 0.025;  // Epsilon
    tmpParam(1) = 1.0;   // True multiscale alpha
    tmpParam(2) = 2.0;   // True multiscale betainv

    std::ofstream output(DATA_PATH + std::string("HomoHomo.txt"), std::ofstream::out | std::ofstream::trunc);
    std::ofstream outputSol(DATA_PATH + std::string("HomoHomoSol.txt"), std::ofstream::out | std::ofstream::trunc);

    // Compute coefficients of the homogenised equation
    std::vector<double> homCoeffs = computeHomCoeffs(tmpParam, (2.0 * M_PI), &V1);
    output << homCoeffs[0] << "\t" << homCoeffs[1] << std::endl;

    // ============= Transform for positiveness ============= //
    VectorXd param = tmpParam;
    param(1) = std::log(param(1));
    param(2) = std::log(param(2));
    // ====================================================== //

    // Initialize structures for the inverse problem
    unsigned long M = 200, nMCMC = 1000;
    double noise = 1e-3;
    double IC = 0.0;
    std::random_device dev;
    std::default_random_engine noiseSeed{1};
    std::default_random_engine proposalSeed{dev()};
    std::default_random_engine acceptanceSeed{dev()};
    std::normal_distribution<double> noiseDistribution(0.0, noise);
    bool IS = true;

    // Generate and perturb observations
    auto x = generateObservations1D(sdeHomo, IC, param, T, N, dev());
    for (auto const &itSol : x)
        outputSol << itSol << "\t";
    outputSol << std::endl;
    outputSol << x[0] << "\t";

    for (unsigned long i = 1; i < x.size(); i++) {
        x[i] += noiseDistribution(noiseSeed);
        // outputSol << x[i] << "\t";
    }
    // outputSol << std::endl;

    // Compute a homogeneous path with the same Brownian for comparison
    VectorXd homoParam(3);
    homoParam(0) = 0.1;
    homoParam(1) = std::log(homCoeffs[0]);
    homoParam(2) = std::log(homCoeffs[1]);
    auto xHom = generateObservations1D(sdeHomo, IC, homoParam, T, N, 0);
    for (auto const &itSol : xHom) {
        outputSol << itSol << "\t";
    }
    outputSol << std::endl;
    outputSol.close();

    // Initial parameter guess
    VectorXd initGuess = VectorXd::Zero(param.size());
    initGuess(0) = param(0);

    system("exec rm -r /net/smana3/vol/vol2/anmc/garegnan/Desktop/Project/SDECode/plots/*");
    // Inverse problem
    std::shared_ptr<Posterior> posterior;
    posterior = std::make_shared<SDEPosterior>(x, T, IC, 1, noise, sdeHomo, param(0), M);
    posterior->computePosterior(param);
    // posterior = std::make_shared<PFPosteriorHom>(x, T, IC, 1, noise, sdeHomo, &V1, param(0), M, IS);
    std::shared_ptr<Proposals> proposal;
    std::vector<double> factors = {1.0, 1.0, 1.0};
    proposal = std::make_shared<Proposals>(1e-1, factors);
    MCMC mcmc(initGuess, proposal, posterior, nMCMC);
    auto sample = mcmc.compute(&proposalSeed, &acceptanceSeed);
    for (auto const &itSample : sample)
        output << itSample.transpose() << std::endl;

    output.close();
    return 0;
}