#include <fstream>
#include <iostream>
#include <iomanip>
#include "generateObservations.hpp"
#include "computeHomogeneous.hpp"
#include <MCMC.hpp>
#include <ParFilLib.hpp>
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
    unsigned int N = 5000;

    VectorXd tmpParam(3);
    tmpParam(0) = 0.05;  // Epsilon
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
    unsigned long M = 40 , nMCMC = 10000;
    double noise = 1e-3;
    double IC = 0.0;
    std::random_device dev;
    std::default_random_engine noiseSeed{1};
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
    VectorXd priorMean = VectorXd::Zero(param.size());
    priorMean(0) = param(0);
    VectorXd initGuess = priorMean;
    std::vector<VectorXd> sample = {};

    sample.push_back(initGuess);
    VectorXd priorStdDev(param.size());
    priorStdDev << 0.0, 1.0, 1.0;
    unsigned int nMC = 4000;
    unsigned int nParam = 20;
    double propStdDev = 2e-2;

    std::vector<double> timeVec(N+1);
    for (unsigned int i = 0; i < N+1; i++) {
        timeVec[i] = T/N * i;
    }

    std::vector<std::vector<VectorXd>> sampleThreads(nParam);
    for (unsigned int init = 0; init < nParam; init++) {
        sampleThreads[init].push_back(initGuess);
    }

    unsigned int L = 20;
    for (unsigned int l = 0; l < L; l++) {
        std::cout << "iteration " << l+1 << " / " << L << std::endl;
        // Compute the modeling error statistics and rescale the observations
        std::vector<std::vector<double>> errors;

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
        std::cout << priorMean.transpose() << std::endl << priorStdDev.transpose() << std::endl;

        std::cout << "Computing modeling error..." << std::endl;
        ModErrAll modErr(sdeHomo, sde, &V1, IC, priorMean, priorStdDev, T, N, x, noise);
        modErr.computePF(nParam, nMC);
        modErr.getModErr(errors);
        std::cout << "Computed modeling error" << std::endl;

        bool plot = true;
        if (plot) {
            plt::named_plot("xe", timeVec, x, "b");
            plt::named_plot("x0", timeVec, xHom, "r");
        }
        std::vector<double> rescaledObs(N + 1);
        for (unsigned int i = 0; i < errors.size(); i++) {
            for (unsigned int j = 0; j < N + 1; j++) {
                rescaledObs[j] = x[j] - errors[i][j];
            }
            if (plot) {
                if (!i) {
                    plt::named_plot("xt", timeVec, rescaledObs, "k");
                } else {
                    plt::plot(timeVec, rescaledObs, "k");
                }
            }
        }
        if (plot) {
            plt::legend();
            plt::show();
        }

        #pragma omp parallel for num_threads(5)
        for (unsigned int modIt = 0; modIt < errors.size(); modIt++) {
            // Inverse problem
            std::vector<double> dummy = {};
            std::shared_ptr<Posterior> posterior;
            std::vector<std::vector<double>> errorsThreads(1, std::vector<double>(N+1));
            for (unsigned int itit = 0; itit < N+1; itit++) {
                errorsThreads[0][itit] = errors[modIt][itit];
            }
            posterior = std::make_shared<PFPosteriorHom>(x, T, IC, 1, noise, sdeHomo, &V1, param(0), M, IS, dummy, &errorsThreads);
            std::shared_ptr<Proposals> proposal;
            std::vector<double> factors = {1.0, 20.0, 1.0};
            /* if (l > 0) {
                factors[1] = priorStdDev(1) / priorStdDev(2);
                propStdDev = priorStdDev(2);
                std::cout << "proposal factors: " << factors[1] << " " << propStdDev << std::endl;
            } */
            proposal = std::make_shared<Proposals>(propStdDev, factors);
            MCMC mcmc(sampleThreads[modIt].back(), proposal, posterior, nMCMC / L);
            sampleThreads[modIt] = mcmc.compute(&proposalSeed, &acceptanceSeed);
        }

        sample.clear();
        for (unsigned int modIt = 0; modIt < errors.size(); modIt++) {
            sample.insert(sample.end(), sampleThreads[modIt].begin(), sampleThreads[modIt].end());
        }
        for (auto const &it : sample) {
            output << it.transpose() << std::endl;
        }
    }

    output.close();
    outputSol.close();
    return 0;
}