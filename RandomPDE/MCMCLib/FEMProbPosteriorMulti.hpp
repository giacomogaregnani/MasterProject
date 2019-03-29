#ifndef FEMProbPosteriorMultiPROB_HPP
#define FEMProbPosteriorMultiPROB_HPP

#include <Eigen/Dense>
#include <OneDimEllipticSolver.hpp>
#include "KarhunenLoeve.hpp"
#include "Proposals.hpp"
#include "Posterior.hpp"

using namespace Eigen;

class FEMProbPosteriorMulti : public Posterior {
private:
    Solver FEMSolver;

    VectorXd originalMesh;

    double p;

    double h;

    int nMC;

    KarhunenLoeve KL;

    VectorXd observations;

    VectorXd xObs;

    double noise;

public:
    FEMProbPosteriorMulti() {};

    FEMProbPosteriorMulti(double (*RHS) (double), double h,
                     std::vector<double>& BC,
                     double p,
                     VectorXd& fieldMean,
                     covKernelFcts  cov,
                     std::vector<double>& KLParam,
                     VectorXd& observations,
                     VectorXd& xObs,
                     double noise);

    double computePosterior(VectorXd& theta);

    VectorXd coeffToField(VectorXd& theta);

    void resetMesh(void);
};

#endif