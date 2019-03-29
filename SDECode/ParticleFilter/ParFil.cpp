#include <iostream>
#include "ParFil.hpp"

ParFil::ParFil(std::vector<double>& y, double T, double IC,
               unsigned int sR, double noise, oneDimSde sde,
               double eps, unsigned long M):
        T(T),
        IC(IC),
        samplingRatio(sR),
        obs(y),
        noise(noise),
        nParticles(M),
        eps(eps)
{
    std::default_random_engine seed{(unsigned int) time(nullptr)};
    Solver = std::make_shared<EM1D>(sde, seed);
    shufflerSeed = std::default_random_engine{(unsigned int) time(nullptr)};
    X.resize(nParticles);
    for (unsigned int i = 0; i < nParticles; i++)
        X[i].resize(obs.size());
    W.resize(nParticles);
    likelihood = 0.0;
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
    auto nObs = static_cast<unsigned int>(std::round(N / samplingRatio));

    Solver->modifyParam(theta);

    unsigned long index = 0;

    for (unsigned long j = 0; j < nObs-1; j++) {

        shuffler = std::discrete_distribution<unsigned int>(W.begin(), W.end());

        for (unsigned long k = 0; k < nParticles; k++)
            X[k][index] = X[shuffler(shufflerSeed)][index];

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

double ParFil::getLikelihood() const
{
    return likelihood;
}

std::vector<double> ParFil::getBestX()
{
    shuffler = std::discrete_distribution<unsigned int>(W.begin(), W.end());
    unsigned int index = shuffler(shufflerSeed);
    return X[index];
}