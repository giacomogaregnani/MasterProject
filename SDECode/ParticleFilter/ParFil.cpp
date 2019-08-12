#include <iostream>
#include <algorithm>
#include "ParFil.hpp"
#include "../matplotlib-cpp-master/matplotlibcpp.h"

namespace plt = matplotlibcpp;

ParFil::ParFil(std::vector<double>& y, double T, double IC, unsigned int sR,
               double noise, oneDimSde &sde, double eps, unsigned long nParticles,
               std::vector<double> timeNoise):
        sde(sde),
        T(T),
        IC(IC),
        samplingRatio(sR),
        obs(y),
        noise(noise),
        nParticles(nParticles),
        eps(eps),
        timeNoise(timeNoise)
{
    std::random_device randomDevice;
    std::default_random_engine EMSeed{randomDevice()};
    Solver = std::make_shared<EM1D>(sde, EMSeed);
    seed = std::default_random_engine{randomDevice()};
    X.resize(nParticles);
    for (unsigned int i = 0; i < nParticles; i++)
        X[i].resize(obs.size());
    XOld.resize(nParticles);
    W.resize(nParticles);
    likelihood = 0.0;
}

void ParFil::compute(VectorXd& theta, std::vector<std::vector<double>>* mod)
{
    bool allModErrFalse = (mod == nullptr);
    bool staticNoise = timeNoise.empty();

    // Keep in memory the Brownian increments
    BM.resize(nParticles);
    for (unsigned int i = 0; i < nParticles; i++)
        BM[i].resize(obs.size()-1);
    BMOld.resize(nParticles);

    // Keep in memory the tree
    tree.resize(nParticles);
    for (unsigned int i = 0; i < nParticles; i++)
        tree[i].resize(obs.size()-1);

    // The first parameter is the multiscale epsilon
    theta(0) = eps;

    for (unsigned int k = 0; k < nParticles; k++) {
        X[k][0] = IC;
        W[k] = 1.0 / nParticles;
    }
    auto N = obs.size();
    double wSum, h = T / (N - 1);
    likelihood = 0;
    auto nObs = (N - 1) / samplingRatio;
    unsigned long index = 0, obsIdx;
    Solver->modifyParam(theta);

    for (unsigned long j = 0; j < nObs; j++) {
        shuffler = std::discrete_distribution<unsigned int>(W.begin(), W.end());
        obsIdx = (j + 1) * samplingRatio;

        for (unsigned long k = 0; k < nParticles; k++) {
            XOld[k] = X[k][index];
            BMOld[k] = BM[k][index];
        }

        for (unsigned long k = 0; k < nParticles; k++) {
            auto shuffleIdx = shuffler(seed);
            tree[k][index] = shuffleIdx;
            X[k][index] = XOld[shuffleIdx];
            BM[k][index] = BMOld[shuffleIdx];
        }

        wSum = 0.0;
        for (unsigned long k = 0; k < nParticles; k++) {
            for (unsigned long i = 0; i < samplingRatio; i++) {
                BM[k][index+i] = std::sqrt(h) * gaussian(seed);
                X[k][index+i+1] = Solver->oneStepGivenNoise(h, X[k][index+i], BM[k][index+i]);
            }
            if (staticNoise) {
                if (allModErrFalse) {
                    W[k] = gaussianDensity(X[k][index + samplingRatio], obs[obsIdx], noise);
                } else {
                    W[k] = 0.0;
                    double temp;
                    for (unsigned long idxmod = 0; idxmod < mod->size(); idxmod++) {
                        temp = obs[obsIdx] - (*mod)[idxmod][obsIdx];
                        W[k] += gaussianDensity(X[k][index + samplingRatio], temp, noise);
                    }
                    W[k] /= mod->size();
                }
            } else {
                double trueStdDev = std::sqrt(noise * noise + timeNoise[obsIdx] * timeNoise[obsIdx]);
                W[k] = gaussianDensity(X[k][index + samplingRatio], obs[obsIdx], trueStdDev);
            }
            wSum += W[k];
        }
        index = index + samplingRatio;

        std::transform(W.begin(), W.end(), W.begin(), [wSum](double& c){return c/wSum;});
        likelihood += std::log(wSum / nParticles);
    }

    // Plot a trajectory
    /* std::vector<double> timeVec(N);
    for (unsigned int i = 0; i < N; i++) {
        timeVec[i] = h*i;
    }
    auto xPlot = sampleX();
    plt::plot(timeVec, obs, "b");
    plt::plot(timeVec, xPlot, "r");
    plt::show(); */
}

double ParFil::importanceSampler(double h, double hObs, double x, VectorXd &theta,
                                 unsigned long obsIdx, unsigned long j, double trueNoise, double correction)
{
    // Notation refers to Golightly, Wilkinson (2011) appendix A.3
    double alphaj = sde.drift(x, theta),
           betaj = sde.diffusion(x, theta),
           deltaj = hObs - j*h;

    if (trueNoise == 0)
        trueNoise = noise;

    // In Golightly, Wilkinson (2011) the diffusion is under square root
    betaj *= betaj;

    double denom = betaj * deltaj + trueNoise * trueNoise;
    double aj = alphaj + betaj * (obs[obsIdx] - correction - (x + alphaj * deltaj)) / denom;
    double bj = betaj - h * betaj * betaj / denom;

    ISmean = x + aj * h;
    ISstddev = std::sqrt(bj * h);

    gaussian.param(std::normal_distribution<double>::param_type(ISmean, ISstddev));

    return gaussian(seed);
}

void ParFil::computeDiffBridge(VectorXd& theta, std::vector<std::vector<double>>* mod,
                               std::vector<double>* weights, bool verbose)
{
    bool allModErrFalse = (mod == nullptr);
    bool staticNoise = timeNoise.empty();
    bool weightsNoiseFalse = (mod == nullptr);

    // The first parameter is the multiscale epsilon
    theta(0) = eps;

    // Initialise
    for (unsigned int k = 0; k < nParticles; k++) {
        X[k][0] = IC;
        W[k] = 1.0 / nParticles;
    }
    auto N = obs.size();
    double h = T / (N-1);

    auto nObs = (N - 1) / samplingRatio;
    unsigned long obsIdx;
    double wSum, hObs = T/nObs, transDensMean, transDensStddev,
           obsDens, transDens, ISDens, temp;
    likelihood = 0.0;
    unsigned long index = 0;

    long maxWeightIdx;
    if (!weightsNoiseFalse)
        maxWeightIdx = std::max_element((*weights).begin(), (*weights).end()) - (*weights).begin();

    for (unsigned long j = 0; j < nObs; j++) {
        shuffler = std::discrete_distribution<unsigned int>(W.begin(), W.end());
        obsIdx = (j + 1) * samplingRatio;

        for (unsigned long k = 0; k < nParticles; k++)
            XOld[k] = X[k][index];

        for (unsigned long k = 0; k < nParticles; k++)
            X[k][index] = XOld[shuffler(seed)];

        wSum = 0.0;
        for (unsigned long k = 0; k < nParticles; k++) {
            transDens = 1.0;
            ISDens = 1.0;
            for (unsigned long i = 0; i < samplingRatio; i++) {
                // Sample from the IS density
                if (staticNoise) {
                    if (allModErrFalse) {
                        temp = importanceSampler(h, hObs, X[k][index + i], theta, obsIdx, i);
                    } else {
                        if (weightsNoiseFalse) {
                            temp = importanceSampler(h, hObs, X[k][index + i], theta, obsIdx, i, 0,
                                                     mod->back()[obsIdx]);
                        } else {
                            temp = importanceSampler(h, hObs, X[k][index + i], theta, obsIdx, i, 0,
                                                     (*mod)[maxWeightIdx][obsIdx]);
                        }
                    }
                } else {
                    double trueStdDev = std::sqrt(noise * noise + timeNoise[obsIdx] * timeNoise[obsIdx]);
                    temp = importanceSampler(h, hObs, X[k][index + i], theta, obsIdx, i, trueStdDev);
                }
                // Evaluate the transition density of the true process
                transDensMean = X[k][index+i] + sde.drift(X[k][index+i], theta) * h;
                transDensStddev = sde.diffusion(X[k][index+i], theta) * std::sqrt(h);
                transDens *= gaussianDensity(temp, transDensMean, transDensStddev);
                // Evaluate the IS density
                ISDens *= gaussianDensity(temp, ISmean, ISstddev);
                X[k][index+1+i] = temp;
            }

            // Evaluate the observation density
            if (staticNoise) {
                if (allModErrFalse) {
                    obsDens = gaussianDensity(X[k][index + samplingRatio], obs[obsIdx], noise);
                } else {
                    if (weightsNoiseFalse) {
                        obsDens = 0.0;
                        double diff;
                        for (unsigned long idxmod = 0; idxmod < mod->size(); idxmod++) {
                            diff = obs[obsIdx] - (*mod)[idxmod][obsIdx];
                            obsDens += gaussianDensity(X[k][index + samplingRatio], diff, noise);
                        }
                        obsDens /= mod->size();
                    } else {
                        obsDens = 0.0;
                        double diff;
                        for (unsigned long idxmod = 0; idxmod < mod->size(); idxmod++) {
                            if ((*weights)[idxmod] > 1e-3) {
                                diff = obs[obsIdx] - (*mod)[idxmod][obsIdx];
                                obsDens += (*weights)[idxmod] * gaussianDensity(X[k][index + samplingRatio], diff, noise);
                            }
                        }
                    }
                }
            } else {
                double trueStdDev = std::sqrt(noise * noise + timeNoise[obsIdx] * timeNoise[obsIdx]);
                obsDens = gaussianDensity(X[k][index + samplingRatio], obs[obsIdx], trueStdDev);
            }

            // Compute the weights
            W[k] = obsDens * transDens / ISDens;
            wSum += W[k];
        }
        index = index + samplingRatio;

        // Normalize and update likelihood
        std::transform(W.begin(), W.end(), W.begin(), [wSum](double& c){return c/wSum;});
        likelihood += std::log(wSum / nParticles);
    }

    if (verbose) {
        double wSumSqd = 0.0;
        for (unsigned long k = 0; k < nParticles; k++) {
            wSumSqd += W[k] * W[k];
        }
        std::cout << "ESS at final time = " << 1.0 / wSumSqd << std::endl;
    }

    // Plot a trajectory
    /* std::vector<double> timeVec(N);
    for (unsigned int i = 0; i < N; i++) {
        timeVec[i] = h*i;
    }
    auto xPlot = sampleX();
    plt::plot(timeVec, obs, "b");
    plt::plot(timeVec, xPlot, "r");
    plt::show(); */
}

double ParFil::getLikelihood() const
{
    return likelihood;
}

std::vector<std::vector<double>> ParFil::getX() const
{
    return X;
}

std::vector<double> ParFil::sampleX()
{
    shuffler = std::discrete_distribution<unsigned int>(W.begin(), W.end());
    unsigned int index = shuffler(seed);
    return X[index];
}

std::vector<std::vector<double>> ParFil::getBM() const
{
    return BM;
}

std::vector<double> ParFil::getW() const
{
    return W;
}

std::vector<std::vector<unsigned int>> ParFil::getTree() const
{
    return tree;
}