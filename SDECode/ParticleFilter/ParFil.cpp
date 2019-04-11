#include <iostream>
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
    std::default_random_engine seed{(unsigned int) time(nullptr)};
    Solver = std::make_shared<EM1D>(sde, seed);
    particleSeed = std::default_random_engine{(unsigned int) time(nullptr)};
    X.resize(nParticles);
    for (unsigned int i = 0; i < nParticles; i++)
        X[i].resize(obs.size());
    W.resize(nParticles);
    likelihood = 0.0;
}


double gaussianDensity(double x, double m, double s)
{
    return 1.0 / (std::sqrt(2.0 * M_PI) * s) * std::exp(-0.5 * (x  - m) * (x - m) / (s * s));
}

void ParFil::compute(VectorXd& theta)
{
    theta(0) = eps;

    // Initialise
    for (unsigned int k = 0; k < nParticles; k++) {
        X[k][0] = IC;
        W[k] = 1.0 / nParticles;
    }

    auto N = obs.size();
    double dist, wSum, h = T/N;
    likelihood = 0;
    auto nObs = static_cast<unsigned int>(std::round((N - 1) / samplingRatio));

    Solver->modifyParam(theta);

    unsigned long index = 0;

    for (unsigned long j = 0; j < nObs; j++) {

        shuffler = std::discrete_distribution<unsigned int>(W.begin(), W.end());

        for (unsigned long k = 0; k < nParticles; k++)
            X[k][index] = X[shuffler(particleSeed)][index];

        wSum = 0.0;
        for (unsigned long k = 0; k < nParticles; k++) {
            for (unsigned long i = 0; i < samplingRatio; i++) {
                X[k][index + 1 + i] = Solver->oneStep(h, X[k][index + i]);
            }
            dist = X[k][index + samplingRatio] - obs[(j + 1) * samplingRatio];
            W[k] = std::exp(-0.5 / (noise * noise) * dist * dist);
            wSum += W[k];
        }
        index = index + samplingRatio;

        for (unsigned long k = 0; k < nParticles; k++) {
            W[k] /= wSum;
        }

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

    double denom = betaj * deltaj + noise * noise;
    double aj = alphaj + betaj * (obs[obsIdx] - (x + alphaj * deltaj)) / denom,
           bj = betaj - betaj / denom * betaj * h;

    ISmean = x + aj * h;
    ISstddev = std::sqrt(bj * h);
    ISgaussian.param(std::normal_distribution<double>::param_type(ISmean, ISstddev));

    return ISgaussian(particleSeed);
}

void ParFil::computeDiffBridge(VectorXd &theta)
{
    theta(0) = eps;

    // Initialise
    for (unsigned int k = 0; k < nParticles; k++) {
        X[k][0] = IC;
        W[k] = 1.0 / nParticles;
    }

    auto N = obs.size();
    auto nObs = static_cast<unsigned int>(std::round((N - 1) / samplingRatio));
    unsigned long obsIdx;
    double wSum, h = T/N, hObs = T/nObs, transDensMean, transDensStddev,
           obsDens, transDens, ISDens, temp;
    likelihood = 0;

    Solver->modifyParam(theta);

    unsigned long index = 0;

    for (unsigned long j = 0; j < nObs; j++) {

        shuffler = std::discrete_distribution<unsigned int>(W.begin(), W.end());
        obsIdx = (j + 1) * samplingRatio;

        for (unsigned long k = 0; k < nParticles; k++)
            X[k][index] = X[shuffler(particleSeed)][index];

        wSum = 0.0;
        for (unsigned long k = 0; k < nParticles; k++) {
            transDens = 1.0, ISDens = 1.0;
            for (unsigned long i = 0; i < samplingRatio; i++) {
                // Sample from the IS density
                temp = importanceSampler(h, hObs, X[k][index + i], theta, obsIdx, i);
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

        for (unsigned long k = 0; k < nParticles; k++) {
            W[k] /= wSum;
        }

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