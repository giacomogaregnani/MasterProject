#include <iostream>
#include "SDEPosterior.hpp"

SDEPosterior::SDEPosterior(std::vector<double>& x, double T, double IC,
                           unsigned int sR, double noise, oneDimSde sde,
                           double eps, unsigned long M):
        T(T),
        IC(IC),
        samplingRatio(sR),
        obs(x),
        noise(noise),
        eps(eps),
        nMC(M)
{
    std::default_random_engine seed{(unsigned int) time(nullptr)};
    Solver = EM1D(sde, seed);
}

double SDEPosterior::computePosterior(VectorXd& theta)
{
    double prior = -0.5 * theta.dot(theta);
    double solution, likelihood, dist;
    auto N = obs.size();
    double h = T/N;
    auto nObs = static_cast<unsigned int>(std::round(N / samplingRatio));

    Solver.modifyParam(theta);

    std::vector<double> likelihoods(nMC);
    for (unsigned long i = 0; i < nMC; i++)
        likelihoods[i] = 0.0;

    #pragma omp parallel for num_threads(6)
    for (unsigned long k = 0; k < nMC; k++) {
        solution = IC;

        for (unsigned long j = 0; j < nObs; j++) {
            for (unsigned long i = 0; i < samplingRatio; i++) {
                solution = Solver.oneStep(h, solution);
            }
            dist = solution - obs[(j+1)*samplingRatio];
            likelihoods[k] += -0.5 / (noise * noise) * dist * dist;
        }
    }

    auto maxLikIt = std::max_element(likelihoods.begin(), likelihoods.end());
    double maxLik = *maxLikIt;
    likelihoods.erase(maxLikIt);

    double sum = 0;
    for (auto it : likelihoods) {
        sum += exp(it - maxLik);
    }
    likelihood = maxLik + std::log(1.0 + sum) - std::log(nMC);

    return prior + likelihood;
}