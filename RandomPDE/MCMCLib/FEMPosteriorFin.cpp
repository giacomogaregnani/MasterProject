#include <iostream>
#include "FEMPosteriorFin.hpp"

FEMPosteriorFin::FEMPosteriorFin(double (*RHS) (double), double h,
                                         std::vector<double>& BC,
                                         VectorXd (*buildField) (VectorXd&, double),
                                         VectorXd& observations,
                                         VectorXd& xObs,
                                         double noise):
        h(h),
        nMC(nMC),
        observations(observations),
        xObs(xObs),
        noise(noise),
        buildField(buildField)
{
    FEMSolver = Solver(RHS, 0.0, 1.0, h, BC[0], BC[1]);
}

double FEMPosteriorFin::computePosterior(VectorXd& theta)
{
    VectorXd field = buildField(theta, h);

    VectorXd uInv = FEMSolver.solve(field);
    VectorXd dist = FEMSolver.evaluate(xObs) - observations;

    double likelihood = -0.5 / (noise * noise) * dist.dot(dist);

    double prior = -0.5 * theta.dot(theta);
    return likelihood + prior;
}

VectorXd FEMPosteriorFin::coeffToField(VectorXd& theta)
{
    VectorXd result = buildField(theta, h);
    return result;
}