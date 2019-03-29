#include "RKEEPosteriorIC.hpp"

RKEEPosteriorIC::RKEEPosteriorIC(double h,
                                 std::vector<VectorXd>& observations,
                                 std::vector<VectorXd>& modErrors,
                                 std::vector<double>& tObs,
                                 double noise,
                                 odeDef ODE,
                                 Butcher tableau):
        ODE(ODE),
        h(h),
        observations(observations),
        tObs(tObs),
        isIdentity(true),
        noise(noise),
        modErrors(modErrors)
{
    RKSolver = RungeKutta(ODE, tableau);
}

RKEEPosteriorIC::RKEEPosteriorIC(double h,
                                 std::vector<VectorXd>& observations,
                                 std::vector<VectorXd>& modErrors,
                                 std::vector<double>& tObs,
                                 MatrixXd& noise,
                                 odeDef ODE,
                                 Butcher tableau):
        ODE(ODE),
        h(h),
        observations(observations),
        tObs(tObs),
        isIdentity(false),
        noiseMatrix(noise),
        modErrors(modErrors)
{
    RKSolver = RungeKutta(ODE, tableau);
    noiseMatrixLLT.compute(noiseMatrix);
}


double RKEEPosteriorIC::computePosterior(VectorXd& theta)
{
    VectorXd RKSolution = theta;
    double prior = -0.5 * theta.dot(theta);

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
    }

    VectorXd dist(modErrors.size());
    std::vector<double> likelihoods(modErrors.size(), 0.0);

    for (unsigned int k = 0; k < modErrors.size(); k++) {
        dist = observations[0] - RKSolution - modErrors[k];
        if (isIdentity) {
            likelihoods[k] += -0.5 / (noise * noise) * dist.dot(dist);
        } else {
            likelihoods[k] += -0.5 * dist.dot(noiseMatrixLLT.solve(dist));
        }
    }

    auto maxLikIt = std::max_element(likelihoods.begin(), likelihoods.end());
    double maxLik = *maxLikIt;
    likelihoods.erase(maxLikIt);

    double sum = 0;
    for (auto it : likelihoods) {
        sum += exp(it - maxLik);
    }
    double likelihood = maxLik + std::log(1.0 + sum) - std::log(modErrors.size());

    return likelihood + prior;
}
