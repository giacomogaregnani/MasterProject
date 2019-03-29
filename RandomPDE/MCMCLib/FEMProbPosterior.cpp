#include <iostream>
#include "FEMProbPosterior.hpp"

FEMProbPosterior::FEMProbPosterior(double (*RHS) (double), double h,
                           std::vector<double>& BC,
                           double p,
                           int nMC,
                           VectorXd& fieldMean,
                           covKernelFcts  cov,
                           std::vector<double>& KLParam,
                           VectorXd& observations,
                           VectorXd& xObs,
                           double noise,
                           bool CNProposal):
        p(p),
        h(h),
        nMC(nMC),
        observations(observations),
        xObs(xObs),
        noise(noise),
        isCN(CNProposal)
{
    FEMSolver = Solver(RHS, 0.0, 1.0, h, BC[0], BC[1]);
    auto N = static_cast<int>(1.0 / h);
    KL = KarhunenLoeve(fieldMean, cov, N, KLParam);
}

double FEMProbPosterior::computePosterior(VectorXd& theta)
{
    std::vector<double> likelihoods(nMC);
    VectorXd mesh = FEMSolver.getMesh();
    long N = mesh.size() - 1;
    VectorXd intPoints = mesh.segment(1, N - 1);
    VectorXd newMesh;

    VectorXd field = KL.KL(theta);

    for (int i = 0; i < field.size(); i++)
        field(i) = std::exp(field(i));

    for (int k = 0; k < nMC; k++) {
        newMesh = mesh;
        newMesh.segment(1, N - 1) += h * h / 2.0 * VectorXd::Random(N - 1);
        FEMSolver.changeMesh(newMesh);
        VectorXd uInv = FEMSolver.solve(field);
        FEMSolver.changeMesh(mesh);
        VectorXd dist = FEMSolver.evaluate(xObs) - observations;
        likelihoods[k] = -0.5 / (noise * noise) * dist.dot(dist);
    }
    FEMSolver.changeMesh(mesh);

    auto maxLikIt = std::max_element(likelihoods.begin(), likelihoods.end());
    double likelihood = std::accumulate(likelihoods.begin(), likelihoods.end(), 0.0) / likelihoods.size();
    /* double maxLik = *maxLikIt;
    likelihoods.erase(maxLikIt);

    double sum = 0;
    for (auto it : likelihoods) {
        sum += exp(it - maxLik);
    }
    double likelihood = maxLik + std::log(1.0 + sum) - std::log(static_cast<double>(nMC)); */

    if (isCN)
        return likelihood;

    double prior = -0.5 * theta.dot(theta);
    return likelihood + prior;
}

VectorXd FEMProbPosterior::coeffToField(VectorXd& theta)
{
    VectorXd result = KL.KL(theta);

    for (int i = 0; i < result.size(); i++)
        result(i) = std::exp(result(i));

    return result;
}