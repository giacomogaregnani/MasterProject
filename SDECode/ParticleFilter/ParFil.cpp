#include <iostream>
#include <algorithm>
#include "ParFil.hpp"
#include "../matplotlib-cpp-master/matplotlibcpp.h"

namespace plt = matplotlibcpp;

ParFil::ParFil(std::vector<double>& y, double T, double IC,
               double noise, double eps, unsigned long nParticles,
               std::shared_ptr<ForwardPF>& forwardSampler):
        T(T),
        IC(IC),
        obs(y),
        noise(noise),
        nParticles(nParticles),
        eps(eps),
        forwardSampler(forwardSampler)
{
    std::random_device randomDevice;
    seed = std::default_random_engine{randomDevice()};
    X.resize(nParticles);
    for (unsigned int i = 0; i < nParticles; i++)
        X[i].resize(obs.size());
    XOld.resize(nParticles);
    W.resize(nParticles);
    likelihood = 0.0;
}

void ParFil::compute(VectorXd& theta, std::vector<std::vector<double>>* mod,
                     std::vector<double>* weights)
{
    bool allModErrFalse = (mod == nullptr);

    if (!allModErrFalse && weights == nullptr) {
        throw std::invalid_argument("Weights have to be provided");
    }

    // The first parameter is the multiscale epsilon
    theta(0) = eps;

    for (unsigned int k = 0; k < nParticles; k++) {
        X[k][0] = IC;
        W[k] = 1.0 / nParticles;
    }
    auto N = obs.size() - 1;
    double wSum;
    likelihood = 0;

    forwardSampler->modifyParam(theta);

    for (unsigned long j = 0; j < N; j++) {
        shuffler = std::discrete_distribution<unsigned int>(W.begin(), W.end());

        for (unsigned long k = 0; k < nParticles; k++)
            XOld[k] = X[k][j];
        for (unsigned long k = 0; k < nParticles; k++) {
            auto shuffleIdx = shuffler(seed);
            X[k][j] = XOld[shuffleIdx];
        }

        wSum = 0.0;
        for (unsigned long k = 0; k < nParticles; k++) {
            X[k][j+1] = forwardSampler->generateSample(X[k][j]);
            if (allModErrFalse) {
                W[k] = gaussianDensity(X[k][j+1], obs[j+1], noise);
            } else {
                W[k] = 0.0;
                double diff;
                for (unsigned long idxmod = 0; idxmod < mod->size(); idxmod++) {
                    if ((*weights)[idxmod] > 1e-3) {
                        diff = obs[j+1] - (*mod)[idxmod][j+1];
                        W[k] += (*weights)[idxmod] * gaussianDensity(X[k][j+1], diff, noise);
                    }
                }
                W[k] /= mod->size();
            }
            wSum += W[k];
        }
        std::transform(W.begin(), W.end(), W.begin(), [wSum](double& c){return c/wSum;});
        likelihood += std::log(wSum / nParticles);
    }
}

void ParFil::computeDiffBridge(VectorXd& theta, std::vector<std::vector<double>>* mod,
                               std::vector<double>* weights, bool verbose)
{
    bool allModErrFalse = (mod == nullptr);

    if (!allModErrFalse && weights == nullptr) {
        throw std::invalid_argument("Weights have to be provided");
    }

    long maxWeightIdx;
    if (!allModErrFalse)
        maxWeightIdx = std::max_element((*weights).begin(), (*weights).end()) - (*weights).begin();

    // The first parameter is the multiscale epsilon
    theta(0) = eps;

    for (unsigned int k = 0; k < nParticles; k++) {
        X[k][0] = IC;
        W[k] = 1.0 / nParticles;
    }
    auto N = obs.size() - 1;
    double wSum;
    likelihood = 0;

    forwardSampler->modifyParam(theta);

    double transDens, ISDens, obsDens;

    for (unsigned long j = 0; j < N; j++) {
        shuffler = std::discrete_distribution<unsigned int>(W.begin(), W.end());

        for (unsigned long k = 0; k < nParticles; k++)
            XOld[k] = X[k][j];
        for (unsigned long k = 0; k < nParticles; k++) {
            auto shuffleIdx = shuffler(seed);
            X[k][j] = XOld[shuffleIdx];
        }

        wSum = 0.0;
        for (unsigned long k = 0; k < nParticles; k++) {
            if (allModErrFalse) {
                X[k][j+1] = forwardSampler->generateSampleIS(X[k][j], obs[j+1], noise);
            } else {
                X[k][j+1] = forwardSampler->generateSampleIS(X[k][j], obs[j+1]-(*mod)[maxWeightIdx][j+1], noise);
            }
            transDens = forwardSampler->evalTransDensity(X[k][j], X[k][j+1]);
            ISDens = forwardSampler->evalISDensity(X[k][j+1]);
            if (allModErrFalse) {
                obsDens = gaussianDensity(X[k][j+1], obs[j+1], noise);
            } else {
                obsDens = 0.0;
                double diff;
                for (unsigned long idxmod = 0; idxmod < mod->size(); idxmod++) {
                    if ((*weights)[idxmod] > 1e-3) {
                        diff = obs[j+1] - (*mod)[idxmod][j+1];
                        obsDens += (*weights)[idxmod] * gaussianDensity(X[k][j+1], diff, noise);
                    }
                }
            }
            W[k] = obsDens * transDens / ISDens;
            wSum += W[k];
        }
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