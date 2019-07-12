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

    double T = 10;
    unsigned int N = 1000;

    VectorXd tmpParam(3);
    tmpParam(0) = 0.2;  // Epsilon
    tmpParam(1) = 1.0;   // True multiscale alpha
    tmpParam(2) = 0.5;   // True multiscale betainv

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
    unsigned long M = 20, nMCMC = 10001;
    double noise = 1e-3;
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

    // Initial parameter guess
    VectorXd priorMean(param.size());
    priorMean << param(0), 0.0, 0.0;
    VectorXd initGuess = priorMean;
    std::vector<VectorXd> sample = {};
    sample.push_back(initGuess);
    VectorXd priorStdDev(param.size());
    priorStdDev << 0.0, 1.0, 1.0;
    unsigned int nMC = 40;
    unsigned int nParam = 10;
    double propStdDev = 2e-2;

    std::vector<double> timeVec(N+1);
    for (unsigned int i = 0; i < N+1; i++) {
        timeVec[i] = T/N * i;
    }

    unsigned int L = 10;
    for (unsigned int l = 0; l < L; l++) {
        // Compute the modeling error statistics and rescale the observations
        std::vector<std::vector<double>> errors;

        std::cout << "Computing modeling error..." << std::endl;
        if (l > 0) {
            priorMean = VectorXd::Zero(3);
            for (auto const &it : sample) {
                priorMean += it / sample.size();
            }
            priorStdDev = VectorXd::Zero(3);
            for (auto const &it : sample) {
                priorStdDev += ((it - priorMean).cwiseProduct(it - priorMean)) / (sample.size() - 1);
            }
            priorStdDev = priorStdDev.array().sqrt();
        }
        ModErrAll modErr(sdeHomo, sde, &V1, IC, priorMean, priorStdDev, T, N, x, noise);
        modErr.computePF(nParam, nMC);
        modErr.getStats(errors);
        std::cout << "Computed modeling error" << std::endl;

        plt::named_plot("xe", timeVec, x, "b");
        plt::named_plot("x0", timeVec, xHom, "r");
        std::vector<double> rescaledObs(N+1);
        for (unsigned int i = 0; i < errors.size()-1; i += 10) {
            for (unsigned int j = 0; j < N+1; j++) {
                rescaledObs[j] = x[j] - errors[i][j];
            }
            if (!i) {
                plt::named_plot("xt", timeVec, rescaledObs, "k");
            } else {
                plt::plot(timeVec, rescaledObs, "k");
            }
        }
        for (unsigned int j = 0; j < N+1; j++) {
            rescaledObs[j] = x[j] - errors.back()[j];
        }
        plt::named_plot("xtm", timeVec, rescaledObs, "g");
        plt::legend();
        plt::show();
        std::cout << priorMean.transpose() << std::endl << priorStdDev.transpose() << std::endl;

        // Inverse problem
        std::vector<double> dummy = {};
        std::shared_ptr<Posterior> posterior;
        // posterior = std::make_shared<PFPosterior>(x, T, IC, 1, noise, sdeHomo, param(0), M, IS); //, stdDevs);
        if (l < 6) {
            posterior = std::make_shared<PFPosteriorHom>(rescaledObs, T, IC, 1, noise, sdeHomo, &V1, param(0), M, IS);
        } else {
            posterior = std::make_shared<PFPosteriorHom>(x, T, IC, 1, noise, sdeHomo, &V1, param(0), M, IS, dummy,
                                                         &errors);
        }
        std::shared_ptr<Proposals> proposal;
        std::vector<double> factors = {1.0, 50.0, 1.0};
        if (l > 0) {
            factors[1] = priorStdDev(1) / priorStdDev(2);
            propStdDev = priorStdDev(2);
            std::cout << "proposal factors: " << factors[1] << " " << propStdDev << std::endl;
        }
        proposal = std::make_shared<Proposals>(propStdDev, factors);
        MCMC mcmc(sample.back(), proposal, posterior, nMCMC/L);
        sample = mcmc.compute(&proposalSeed, &acceptanceSeed);

        for (auto const &it : sample) {
            output << it.transpose() << std::endl;
        }
    }

    output.close();
    outputSol.close();
    return 0;
}