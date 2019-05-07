#include "RKAddNoisePosteriorIC.hpp"

RKAddPosteriorIC::RKAddPosteriorIC(double h,
                                     std::vector<VectorXd>& observations,
                                     std::vector<double>& tObs,
                                     double noise,
                                     odeDef ODE,
                                     Butcher tableau,
                                     int nMC, double p):
        ODE(ODE),
        h(h),
        observations(observations),
        tObs(tObs),
        noise(noise),
        nMC(nMC),
        p(p)
{
    RKSolvers.resize(nMC);
    generators.resize(nMC);
    for (int i = 0; i < nMC; i++) {
        generators[i] = std::default_random_engine{(unsigned int) time(nullptr) + i};
        RKSolvers[i] = RungeKuttaAddNoise(&generators[i], ODE, tableau, h, p);
    }
}

double RKAddPosteriorIC::computePosterior(VectorXd& theta)
{
    double prior = -0.5 * theta.dot(theta);

    // Likelihood computation
    std::vector<double> likelihoods(nMC, 0.0);

    int k = 0;
    for (k = 0; k < nMC; k++) {

        VectorXd RKSolution = theta;
        VectorXd dist = RKSolution;

        for (unsigned int i = 0; i < tObs.size(); i++) {
            double dT;
            if (i == 0) {
                dT = tObs[i];
            } else {
                dT = tObs[i] - tObs[i - 1];
            }
            auto nSteps = static_cast<unsigned long>(std::round(dT / h));

            for (unsigned long j = 0; j < nSteps; j++) {
                RKSolution = RKSolvers[k].oneStep(RKSolution, ODE.refParam);
            }

            dist = RKSolution - observations[i];
            likelihoods[k] += -0.5 / (noise * noise) * dist.dot(dist);
        }
    }

    auto maxLikIt = std::max_element(likelihoods.begin(), likelihoods.end());
    double maxLik = *maxLikIt;
    likelihoods.erase(maxLikIt);

    double sum = 0;
    for (auto it : likelihoods) {
        sum += std::exp(it - maxLik);
    }
    double likelihood = maxLik + std::log(1.0 + sum) - std::log(nMC);

    return likelihood + prior;
}
