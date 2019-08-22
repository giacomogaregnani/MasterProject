#include <fstream>
#include <iostream>
#include <iomanip>
#include "generateObservations.hpp"
#include "computeHomogeneous.hpp"
#include <MCMC.hpp>
#include <ParFilLib.hpp>

/* #include "../matplotlib-cpp-master/matplotlibcpp.h"
namespace plt = matplotlibcpp; */

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

    double T = 40;
    unsigned int N = 4000;

    VectorXd tmpParam(3);
    tmpParam(0) = 0.1;  // Epsilon
    tmpParam(1) = 1.0;   // True multiscale alpha
    tmpParam(2) = 0.5;   // True multiscale betainv

    // Compute coefficients of the homogenised equation
    std::vector<double> homCoeffs = computeHomCoeffs(tmpParam, (2.0 * M_PI), &V1);

    // ============= Transform for positiveness ============= //
    VectorXd param = tmpParam;
    param(1) = std::log(param(1));
    param(2) = std::log(param(2));
    // ====================================================== //

    // Initialize structures for the inverse problem
    double noise = 1e-2;
    double IC = 0.0;
    std::random_device dev;
    std::default_random_engine noiseSeed{dev()};
    std::normal_distribution<double> noiseDistribution(0.0, noise);

    // Generate and perturb observations
    auto x = generateObservations1D(sde, IC, param, T, N, 0);
    for (unsigned long i = 1; i < x.size(); i++) {
        x[i] += noiseDistribution(noiseSeed);
    }

    // Compute a homogeneous path with the same Brownian for comparison
    VectorXd homoParam(3);
    homoParam(0) = 0.1;
    homoParam(1) = std::log(homCoeffs[0]);
    homoParam(2) = std::log(homCoeffs[1]);
    auto xHom = generateObservations1D(sdeHomo, IC, homoParam, T, N, 0);

    // Initial parameter guess
    VectorXd initGuess = param; // VectorXd::Zero(3);
    initGuess(0) = param(0);
    std::vector<double> modErr(N+1), xPF(N+1), modErrPF(N+1);
    unsigned int nMCMod = 100;
    double h = T / N;
    std::default_random_engine seed(dev());

    std::shared_ptr<ForwardPFModErr> FwdModErr;
    FwdModErr = std::make_shared<ForwardPFModErr>(sde, sdeHomo, h, seed, &V1);
    ParFilMod PFMod(x, T, IC, noise, param(0), nMCMod, FwdModErr);

    std::vector<double> timeVec(N + 1);
    for (unsigned int i = 0; i < N + 1; i++) {
        timeVec[i] = h * i;
        modErr[i] = x[i] - xHom[i];
    }

    unsigned long nParam = 10;
    for (unsigned int j = 0; j < nParam; j++) {
        VectorXd r = VectorXd::Random(3);
        VectorXd theta = initGuess + 0.01 * r;
        PFMod.compute(theta, true);

        auto X = PFMod.sampleX();

        for (unsigned int i = 0; i < N+1; i++) {
            xPF[i] = X[i](0);
            modErrPF[i] = X[i](1);
        }
    }
    /* plt::subplot(2, 1, 1);
    plt::named_plot("x", timeVec, x, "b");
    plt::legend();
    plt::subplot(2, 1, 2);
    plt::named_plot("mod", timeVec, modErr, "b");
    plt::legend();
    plt::show(); */

    return 0;
}