#include <iostream>
#include "FEMProbPosteriorFin.hpp"

FEMProbPosteriorFin::FEMProbPosteriorFin(double (*RHS) (double), double h,
                                             std::vector<double>& BC,
                                             double p,
                                             VectorXd (*buildField) (VectorXd&, double),
                                             VectorXd& observations,
                                             VectorXd& xObs,
                                             double noise):
        p(p),
        h(h),
        nMC(nMC),
        observations(observations),
        xObs(xObs),
        noise(noise),
        buildField(buildField)
{
    FEMSolver = Solver(RHS, 0.0, 1.0, h, BC[0], BC[1]);
    originalMesh = FEMSolver.getMesh();
    auto N = originalMesh.size() - 1;
    VectorXd newMesh = originalMesh;
    newMesh.segment(1, N - 1) += 0.5 * std::pow(h, p) * VectorXd::Random(N - 1);
    FEMSolver.changeMesh(newMesh);
}

double FEMProbPosteriorFin::computePosterior(VectorXd& theta)
{
    VectorXd field = buildField(theta, h);

    VectorXd uInv = FEMSolver.solve(field);
    VectorXd dist = FEMSolver.evaluate(xObs) - observations;

    double likelihood = -0.5 / (noise * noise) * dist.dot(dist);

    double prior = -0.5 * theta.dot(theta);
    return likelihood + prior;
}

VectorXd FEMProbPosteriorFin::coeffToField(VectorXd& theta)
{
    VectorXd result = buildField(theta, h);
    return result;
}

void FEMProbPosteriorFin::resetMesh()
{
    auto N = originalMesh.size() - 1;
    VectorXd newMesh = originalMesh;
    newMesh.segment(1, N - 1) += 0.5 * std::pow(h, p) * VectorXd::Random(N - 1);
    FEMSolver.changeMesh(newMesh);
}