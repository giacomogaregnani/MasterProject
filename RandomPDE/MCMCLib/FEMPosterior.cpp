#include "FEMPosterior.hpp"

FEMPosterior::FEMPosterior(double (*RHS) (double), double h,
                           std::vector<double>& BC,
                           VectorXd& fieldMean,
                           covKernelFcts  cov,
                           std::vector<double>& KLParam,
                           VectorXd& observations,
                           VectorXd& xObs,
                           double noise,
                           bool CNProposal):
    observations(observations),
    xObs(xObs),
    noise(noise),
    isCN(CNProposal)
{
    FEMSolver = Solver(RHS, 0.0, 1.0, h, BC[0], BC[1]);
    int N = static_cast<int>(1.0 / h);
    KL = KarhunenLoeve(fieldMean, cov, N+1, KLParam);
}

double FEMPosterior::computePosterior(VectorXd& theta)
{
    VectorXd field = KL.KL(theta);

    for (int i = 0; i < field.size(); i++)
        field(i) = std::exp(field(i));

    VectorXd uInv = FEMSolver.solve(field);
    VectorXd dist = FEMSolver.evaluate(xObs) - observations;

    double likelihood = -0.5 / (noise * noise) * dist.dot(dist);

    if (isCN)
        return likelihood;

    double prior = -0.5 * theta.dot(theta);
    return likelihood + prior;
}

VectorXd FEMPosterior::coeffToField(VectorXd& theta)
{
    VectorXd result = KL.KL(theta);

    for (int i = 0; i < result.size(); i++)
        result(i) = std::exp(result(i));
    return result;
}