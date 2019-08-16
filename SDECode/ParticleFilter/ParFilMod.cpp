#include <iostream>
#include <algorithm>
#include "ParFilMod.hpp"

ParFilMod::ParFilMod(std::vector<double>& y, double T, double IC,
                       double noise, double eps, unsigned long nParticles,
                       std::shared_ptr<ForwardPFModErr>& forwardSampler):
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
    for (unsigned int i = 0; i < nParticles; i++) {
        X[i].resize(obs.size());
    }
    XOld.resize(nParticles);
    W.resize(nParticles);
    likelihood = 0.0;
}

void ParFilMod::compute(VectorXd& theta, bool verbose)
{
    // The first parameter is the multiscale epsilon
    theta(0) = eps;

    for (unsigned int k = 0; k < nParticles; k++) {
        X[k][0](0) = IC;
        X[k][0](1) = 0;
        W[k] = 1.0 / nParticles;
    }
    auto N = obs.size() - 1;
    double wSum, wSumSqd;
    likelihood = 0;

    forwardSampler->modifyParam(theta);

    for (unsigned long j = 0; j < N; j++) {
        shuffler = std::discrete_distribution<unsigned int>(W.begin(), W.end());
        for (unsigned long k = 0; k < nParticles; k++) {
            XOld[k] = X[k][j];
        }
        for (unsigned long k = 0; k < nParticles; k++) {
            auto shuffleIdx = shuffler(seed);
            X[k][j] = XOld[shuffleIdx];
        }

        wSum = 0.0;
        for (unsigned long k = 0; k < nParticles; k++) {
            X[k][j+1] = forwardSampler->generateSample(X[k][j]);
            W[k] = gaussianDensity(X[k][j+1](0), obs[j+1], noise);
            wSum += W[k];
        }

        std::transform(W.begin(), W.end(), W.begin(), [wSum](double& c){return c/wSum;});
        likelihood += std::log(wSum / nParticles);
    }

    if (verbose) {
        wSumSqd = 0.0;
        for (unsigned long k = 0; k < nParticles; k++) {
            wSumSqd += W[k] * W[k];
        }
        std::cout << "ESS at final time = " << 1.0 / wSumSqd << std::endl;
    }
}

void ParFilMod::computeIS(VectorXd& theta, bool verbose)
{
    // The first parameter is the multiscale epsilon
    theta(0) = eps;

    for (unsigned int k = 0; k < nParticles; k++) {
        X[k][0](0) = IC;
        X[k][0](1) = 0.0;
        W[k] = 1.0 / nParticles;
    }
    auto N = obs.size() - 1;
    double wSum, wSumSqd, transDens, ISDens, obsDens;
    likelihood = 0;

    forwardSampler->modifyParam(theta);

    for (unsigned long j = 0; j < N; j++) {
        shuffler = std::discrete_distribution<unsigned int>(W.begin(), W.end());
        for (unsigned long k = 0; k < nParticles; k++) {
            XOld[k] = X[k][j];
        }
        for (unsigned long k = 0; k < nParticles; k++) {
            auto shuffleIdx = shuffler(seed);
            X[k][j] = XOld[shuffleIdx];
        }

        wSum = 0.0;
        for (unsigned long k = 0; k < nParticles; k++) {
            X[k][j+1] = forwardSampler->generateSampleIS(X[k][j], obs[j+1], noise);
            obsDens = gaussianDensity(X[k][j+1](0), obs[j+1], noise);
            transDens = forwardSampler->evalTransDensity(X[k][j], X[k][j+1]);
            ISDens = forwardSampler->evalISDensity(X[k][j+1]);
            W[k] = obsDens * transDens / ISDens;
            wSum += W[k];
        }

        std::transform(W.begin(), W.end(), W.begin(), [wSum](double& c){return c/wSum;});
        likelihood += std::log(wSum / nParticles);
    }

    if (verbose) {
        wSumSqd = 0.0;
        for (unsigned long k = 0; k < nParticles; k++) {
            wSumSqd += W[k] * W[k];
        }
        std::cout << "ESS at final time = " << 1.0 / wSumSqd << std::endl;
    }
}

double ParFilMod::getLikelihood() const
{
    return likelihood;
}

std::vector<std::vector<Vector2d>> ParFilMod::getX() const
{
    return X;
}

std::vector<Vector2d> ParFilMod::sampleX()
{
    shuffler = std::discrete_distribution<unsigned int>(W.begin(), W.end());
    unsigned int index = shuffler(seed);
    return X[index];
}

std::vector<double> ParFilMod::getW() const
{
    return W;
}
