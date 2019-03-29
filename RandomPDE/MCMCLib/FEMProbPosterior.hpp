#ifndef FEMProbPosteriorPROB_HPP
#define FEMProbPosteriorPROB_HPP

#include <Eigen/Dense>
#include <OneDimEllipticSolver.hpp>
#include "KarhunenLoeve.hpp"
#include "Proposals.hpp"
#include "Posterior.hpp"

using namespace Eigen;

class FEMProbPosterior : public Posterior {
private:
    Solver FEMSolver;

    double p;

    double h;

    int nMC;

    KarhunenLoeve KL;

    VectorXd observations;

    VectorXd xObs;

    double noise;

    bool isCN;

public:
    FEMProbPosterior() {};

    FEMProbPosterior(double (*RHS) (double), double h,
                     std::vector<double>& BC,
                     double p,
                     int nMC,
                     VectorXd& fieldMean,
                     covKernelFcts  cov,
                     std::vector<double>& KLParam,
                     VectorXd& observations,
                     VectorXd& xObs,
                     double noise,
                     bool CNProposal);

    double computePosterior(VectorXd& theta);

    VectorXd coeffToField(VectorXd& theta);
};

#endif