#include "RKPosteriorIC.hpp"

RKPosteriorIC::RKPosteriorIC(double h,
                         std::vector<VectorXd>& observations,
                         std::vector<double>& tObs,
                         double noise,
                         odeDef ODE,
                         Butcher tableau):
        ODE(ODE),
        h(h),
        observations(observations),
        tObs(tObs),
        isIdentity(true),
        noise(noise)
{
    RKSolver = RungeKutta(ODE, tableau);
}

RKPosteriorIC::RKPosteriorIC(double h,
                             std::vector<VectorXd>& observations,
                             std::vector<double>& tObs,
                             MatrixXd noise,
                             odeDef ODE,
                             Butcher tableau):
        ODE(ODE),
        h(h),
        observations(observations),
        tObs(tObs),
        isIdentity(false),
        noiseMatrix(noise)
{
    RKSolver = RungeKutta(ODE, tableau);
    noiseMatrixLLT.compute(noiseMatrix);
}


double RKPosteriorIC::computePosterior(VectorXd& theta)
{
    VectorXd RKSolution = theta;

    double prior = -0.5 * theta.dot(theta);

    double likelihood = 0;
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
            RKSolution = RKSolver.oneStep(h, RKSolution, ODE.refParam);
        }

        dist = RKSolution - observations[i];
        if (isIdentity) {
            likelihood += -0.5 / (noise * noise) * dist.dot(dist);
        } else {
            likelihood += -0.5 * dist.dot(noiseMatrixLLT.solve(dist));
        }
    }

    return likelihood + prior;
}
