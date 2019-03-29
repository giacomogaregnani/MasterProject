#include "RKPosterior.hpp"

RKPosterior::RKPosterior(double h,
                         double T,
                         VectorXd& initialCondition,
                         std::vector<VectorXd>& observations,
                         std::vector<double>& tObs,
                         double noise,
                         odeDef ODE,
                         Butcher tableau):
        h(h),
        T(T),
        IC(initialCondition),
        observations(observations),
        tObs(tObs),
        noise(noise)
{
    RKSolver = RungeKutta(ODE, tableau);
}

double RKPosterior::computePosterior(VectorXd& theta)
{
    // Prior computation
    double prior = -0.5 * theta.dot(theta);

    VectorXd RKSolution = IC;
    double likelihood = 0;
    VectorXd dist = RKSolution;

    for (unsigned int i = 0; i < tObs.size(); i++) {

        double dT;
        if (i == 0) {
            dT = tObs[i];
        } else {
            dT = tObs[i] - tObs[i - 1];
        }
        unsigned long nSteps = static_cast<unsigned long>(std::round(dT / h));

        for (unsigned long j = 0; j < nSteps; j++) {
            RKSolution = RKSolver.oneStep(h, RKSolution, theta);
        }

        dist = RKSolution - observations[i];
        likelihood += -0.5 / (noise * noise) * dist.dot(dist);
    }

    return prior + likelihood;
}