#ifndef FEMPOSTERIOR_HPP
#define FEMPOSTERIOR_HPP

#include <Eigen/Dense>
#include <OneDimEllipticSolver.hpp>
#include "KarhunenLoeve.hpp"
#include "Proposals.hpp"
#include "Posterior.hpp"

using namespace Eigen;

class FEMPosterior : public Posterior {
private:
    Solver FEMSolver;

    KarhunenLoeve KL;

    VectorXd observations;

    VectorXd xObs;

    double noise;

    bool isCN;

public:
    FEMPosterior() = default;

    FEMPosterior(double (*RHS) (double), double h,
                 std::vector<double>& BC,
                 VectorXd& fieldMean,
                 covKernelFcts  cov,
                 std::vector<double>& KLParam,
                 VectorXd& observations,
                 VectorXd& xObs,
                 double noise,
                 bool CNProposal = false);

    double computePosterior(VectorXd& theta);

    VectorXd coeffToField(VectorXd& theta);
};

#endif