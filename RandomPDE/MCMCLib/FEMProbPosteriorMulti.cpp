#include <iostream>
#include "FEMProbPosteriorMulti.hpp"

FEMProbPosteriorMulti::FEMProbPosteriorMulti(double (*RHS) (double), double h,
                                   std::vector<double>& BC,
                                   double p,
                                   VectorXd& fieldMean,
                                   covKernelFcts  cov,
                                   std::vector<double>& KLParam,
                                   VectorXd& observations,
                                   VectorXd& xObs,
                                   double noise):
        p(p),
        h(h),
        nMC(nMC),
        observations(observations),
        xObs(xObs),
        noise(noise)
{
    FEMSolver = Solver(RHS, 0.0, 1.0, h, BC[0], BC[1]);
    originalMesh = FEMSolver.getMesh();
    auto N = originalMesh.size() - 1;
    KL = KarhunenLoeve(fieldMean, cov, N+1, KLParam);
    VectorXd newMesh = originalMesh;
    newMesh.segment(1, N - 1) += 0.5 * std::pow(h, p) * VectorXd::Random(N - 1);
    FEMSolver.changeMesh(newMesh);
}

double FEMProbPosteriorMulti::computePosterior(VectorXd& theta)
{
    VectorXd field = KL.KL(theta);

    for (int i = 0; i < field.size(); i++)
        field(i) = std::exp(field(i));

    VectorXd uInv = FEMSolver.solve(field);
    VectorXd dist = FEMSolver.evaluate(xObs) - observations;

    double likelihood = -0.5 / (noise * noise) * dist.dot(dist);

    double prior = -0.5 * theta.dot(theta);
    return likelihood + prior;
}

VectorXd FEMProbPosteriorMulti::coeffToField(VectorXd& theta)
{
    VectorXd result = KL.KL(theta);

    for (int i = 0; i < result.size(); i++)
        result(i) = std::exp(result(i));

    return result;
}

void FEMProbPosteriorMulti::resetMesh()
{
    auto N = originalMesh.size() - 1;
    VectorXd newMesh = originalMesh;
    newMesh.segment(1, N - 1) += 0.5 * std::pow(h, p) * VectorXd::Random(N - 1);
    FEMSolver.changeMesh(newMesh);
}