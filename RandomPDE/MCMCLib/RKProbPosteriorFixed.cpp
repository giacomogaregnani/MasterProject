#include "RKProbPosteriorFixed.hpp"

RKProbPosteriorFixed::RKProbPosteriorFixed(double h,
                         double T,
                         VectorXd& initialCondition,
                         std::vector<VectorXd>& observations,
                         std::vector<double>& tObs,
                         double noise,
                         odeDef ODE,
                         Butcher tableau,
                         double p, int init):
        h(h),
        T(T),
        IC(initialCondition),
        observations(observations),
        tObs(tObs),
        noise(noise)
{
    RKSolver = RungeKutta(ODE, tableau);
    std::uniform_real_distribution<double> unifDistribution(h - std::pow(h, p), h + std::pow(h, p));
    std::default_random_engine generator{(unsigned long) time(nullptr) + init};
    auto N = static_cast<unsigned long>(std::round(T / h));
    timeSteps = std::vector<double>{};

    for (unsigned long i = 0; i < N; i++) {
        timeSteps.push_back(unifDistribution(generator));
    }
}

double RKProbPosteriorFixed::computePosterior(VectorXd& theta)
{
    // Prior computation
    double prior = -0.5 * theta.dot(theta);

    VectorXd RKSolution = IC;
    double likelihood = 0;
    VectorXd dist = RKSolution;

    unsigned int k = 0;

    for (unsigned int i = 0; i < tObs.size(); i++) {

        double dT;
        if (i == 0) {
            dT = tObs[i];
        } else {
            dT = tObs[i] - tObs[i - 1];
        }
        auto nSteps = static_cast<unsigned long>(std::round(dT / h));

        for (unsigned long j = 0; j < nSteps; j++) {
            RKSolution = RKSolver.oneStep(timeSteps[k++], RKSolution, theta);
        }

        dist = RKSolution - observations[i];
        likelihood += -0.5 / (noise * noise) * dist.dot(dist);
    }

    return prior + likelihood;
}