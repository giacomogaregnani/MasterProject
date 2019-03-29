#ifndef FEMProbPosteriorFinPROB_HPP
#define FEMProbPosteriorFinPROB_HPP

#include <Eigen/Dense>
#include <OneDimEllipticSolver.hpp>
#include "KarhunenLoeve.hpp"
#include "Proposals.hpp"
#include "Posterior.hpp"

using namespace Eigen;

class FEMProbPosteriorFin : public Posterior {
private:
    Solver FEMSolver;

    VectorXd originalMesh;

    double p;

    double h;

    int nMC;

    VectorXd observations;

    VectorXd xObs;

    double noise;

    VectorXd (*buildField) (VectorXd&, double);

public:
    FEMProbPosteriorFin() {};

    FEMProbPosteriorFin(double (*RHS) (double), double h,
                     std::vector<double>& BC,
                     double p,
                     VectorXd (*buildField) (VectorXd&, double),
                     VectorXd& observations,
                     VectorXd& xObs,
                     double noise);

    double computePosterior(VectorXd& theta);

    VectorXd coeffToField(VectorXd& theta);

    void resetMesh(void);
};

#endif