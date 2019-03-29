#include <iomanip>
#include "RKProbPosterior.hpp"

RKProbPosterior::RKProbPosterior(double h,
                                 double T,
                                 VectorXd& initialCondition,
                                 std::vector<VectorXd>& observations,
                                 std::vector<double>& tObs,
                                 double noise,
                                 odeDef ODE,
                                 Butcher tableau,
                                 int nMC, double p):
        h(h),
        T(T),
        IC(initialCondition),
        observations(observations),
        tObs(tObs),
        noise(noise),
        nMC(nMC),
        p(p)
{
    generators.resize(nMC);
    RKSolvers.resize(nMC);
    for (int i = 0; i < nMC; i++) {
        generators[i] = std::default_random_engine{(unsigned int) time(NULL) + i};
        RKSolvers[i] = RungeKuttaAddNoise(&generators[i], ODE, tableau, h, p);
    }
}

double RKProbPosterior::computePosterior(VectorXd& theta)
{
    // Prior computation
    double prior = -0.5 * theta.dot(theta);

    // Likelihood computation
    std::vector<double> likelihoods(nMC, 0.0);

    // #pragma omp parallel for num_threads(10)
    for (int k = 0; k < nMC; k++) {
        VectorXd RKSolution = IC;
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
                RKSolution = RKSolvers[k].oneStep(RKSolution, theta);
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
        sum += exp(it - maxLik);
    }
    double likelihood = maxLik + std::log(1.0 + sum) - std::log(nMC);

    /* std::vector<VectorXd> RKSolutions(nMC, IC);

    double likelihood = 0.0;

    for (unsigned int i = 0; i < tObs.size(); i++) {

        double dT;
        if (i == 0) {
            dT = tObs[i];
        } else {
            dT = tObs[i] - tObs[i - 1];
        }
        auto nSteps = static_cast<unsigned long>(std::round(dT / h));

        VectorXd meanSolution = VectorXd::Zero(IC.size());

        for (int k = 0; k < nMC; k++) {
            for (unsigned long j = 0; j < nSteps; j++) {
                RKSolutions[k] = RKSolvers[k].oneStep(RKSolutions[k], theta);
            }
            meanSolution += RKSolutions[k];
        }
        meanSolution = meanSolution / nMC;

        VectorXd dist = meanSolution - observations[i];
        likelihood += -0.5 / (noise * noise) * dist.dot(dist);
    } */

    return prior + likelihood;
}

