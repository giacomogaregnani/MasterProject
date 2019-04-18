#include <iostream>
#include <algorithm>
#include "ParFil.hpp"

ParFil::ParFil(std::vector<double>& y, double T, double IC, unsigned int sR,
               double noise, oneDimSde &sde, double eps, unsigned long nParticles):
        sde(sde),
        T(T),
        IC(IC),
        samplingRatio(sR),
        obs(y),
        noise(noise),
        nParticles(nParticles),
        eps(eps)
{
    std::random_device randomDevice;
    std::default_random_engine seed{randomDevice()};
    Solver = std::make_shared<EM1D>(sde, seed);
    particleSeed = std::default_random_engine{randomDevice()};
    X.resize(nParticles);
    for (unsigned int i = 0; i < nParticles; i++)
        X[i].resize(obs.size());
    W.resize(nParticles);
    likelihood = 0.0;
}


double gaussianDensity(double x, double mu, double sigma)
{
    return 1.0 / (std::sqrt(2.0 * M_PI) * sigma) * std::exp(-0.5 * (x  - mu) * (x - mu) / (sigma * sigma));
}

void ParFil::compute(VectorXd& theta)
{
    // The first parameter is the multiscale epsilon
    theta(0) = eps;

    // Initialize
    for (unsigned int k = 0; k < nParticles; k++) {
        X[k][0] = IC;
        W[k] = 1.0 / nParticles;
    }
    auto N = obs.size();
    double wSum, h = T/N;
    likelihood = 0;
    auto nObs = (N - 1) / samplingRatio;
    unsigned long index = 0, obsIdx;
    Solver->modifyParam(theta);

    for (unsigned long j = 0; j < nObs; j++) {
        // Sample the particles
        shuffler = std::discrete_distribution<unsigned int>(W.begin(), W.end());
        obsIdx = (j + 1) * samplingRatio;
        for (unsigned long k = 0; k < nParticles; k++)
            X[k][index] = X[shuffler(particleSeed)][index];

        // Innovate and compute the weights
        wSum = 0.0;
        for (unsigned long k = 0; k < nParticles; k++) {
            for (unsigned long i = 0; i < samplingRatio; i++)
                X[k][index+1+i] = Solver->oneStep(h, X[k][index+i]);
            W[k] = gaussianDensity(X[k][index+samplingRatio], obs[obsIdx], noise);
            wSum += W[k];
        }
        index = index + samplingRatio;

        // Normalize and update likelihood
        std::transform(W.begin(), W.end(), W.begin(), [wSum](double& c){return c/wSum;});
        likelihood += std::log(wSum / nParticles);
    }
}

double ParFil::importanceSampler(double h, double hObs, double x, VectorXd &theta,
                                 unsigned long obsIdx, unsigned long j)
{
    // Notation refers to Golightly, Wilkinson (2011) appendix A.3
    double alphaj = sde.drift(x, theta),
           betaj = sde.diffusion(x, theta),
           deltaj = hObs - j*h;

    // In Golightly, Wilkinson (2011) the diffusion is under square root
    betaj *= betaj;

    double denom = betaj * deltaj + noise * noise;
    double aj = alphaj + betaj * (obs[obsIdx] - (x + alphaj * deltaj)) / denom;
    double bj = betaj - h * betaj * betaj / denom;

    ISmean = x + aj * h;
    ISstddev = std::sqrt(bj * h);
    ISgaussian.param(std::normal_distribution<double>::param_type(ISmean, ISstddev));

    return ISgaussian(particleSeed);
}

void ParFil::computeDiffBridge(VectorXd& theta)
{
    // The first parameter is the multiscale epsilon
    theta(0) = eps;

    // Initialise
    for (unsigned int k = 0; k < nParticles; k++) {
        X[k][0] = IC;
        W[k] = 1.0 / nParticles;
    }
    auto N = obs.size();
    double h = T / N;
    // ========== THIS IS JUST A TEST ============ //
    // The last parameter is the sampling dt
    if (theta.size() > 3) {
        double delta = std::exp(theta(theta.size()-1));
        auto ratio = static_cast<unsigned int>(delta / h);
        if (ratio < 1) {
            samplingRatio = 1;
        } else if (ratio > N) {
            samplingRatio = N-1;
        } else {
            samplingRatio = ratio;
        }
        theta(3) = std::log(samplingRatio * h);
    }
    // =========================================== //
    auto nObs = (N - 1) / samplingRatio;
    unsigned long obsIdx;
    double wSum, hObs = T/nObs, transDensMean, transDensStddev,
           obsDens, transDens, ISDens, temp;
    likelihood = 0;
    unsigned long index = 0;

    for (unsigned long j = 0; j < nObs; j++) {
        shuffler = std::discrete_distribution<unsigned int>(W.begin(), W.end());
        obsIdx = (j + 1) * samplingRatio;

        for (unsigned long k = 0; k < nParticles; k++)
            X[k][index] = X[shuffler(particleSeed)][index];

        wSum = 0.0;
        for (unsigned long k = 0; k < nParticles; k++) {
            transDens = 1.0;
            ISDens = 1.0;
            for (unsigned long i = 0; i < samplingRatio; i++) {
                // Sample from the IS density
                temp = importanceSampler(h, hObs, X[k][index+i], theta, obsIdx, i);
                // Evaluate the transition density of the true process
                transDensMean = X[k][index+i] + sde.drift(X[k][index+i], theta) * h;
                transDensStddev = sde.diffusion(X[k][index+i], theta) * std::sqrt(h);
                transDens *= gaussianDensity(temp, transDensMean, transDensStddev);
                // Evaluate the IS density
                ISDens *= gaussianDensity(temp, ISmean, ISstddev);
                X[k][index+1+i] = temp;
            }
            // Evaluate the observation density
            obsDens = gaussianDensity(X[k][index + samplingRatio], obs[obsIdx], noise);
            // Compute the weights
            W[k] = obsDens * transDens / ISDens;
            wSum += W[k];
        }
        index = index + samplingRatio;

        // Normalize and update likelihood
        std::transform(W.begin(), W.end(), W.begin(), [wSum](double& c){return c/wSum;});
        likelihood += std::log(wSum / nParticles);
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

std::vector<double> ParFil::getBestX()
{
    shuffler = std::discrete_distribution<unsigned int>(W.begin(), W.end());
    unsigned int index = shuffler(particleSeed);
    return X[index];
}