#ifndef FEMPosteriorFinPROB_HPP
#define FEMPosteriorFinPROB_HPP

#include <Eigen/Dense>
#include <OneDimEllipticSolver.hpp>
#include "KarhunenLoeve.hpp"
#include "Proposals.hpp"
#include "Posterior.hpp"

using namespace Eigen;

class FEMPosteriorFin : public Posterior {
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
    FEMPosteriorFin() {};

    FEMPosteriorFin(double (*RHS) (double), double h,
                        std::vector<double>& BC,
                        VectorXd (*buildField) (VectorXd&, double),
                        VectorXd& observations,
                        VectorXd& xObs,
                        double noise);

    double computePosterior(VectorXd& theta);

    VectorXd coeffToField(VectorXd& theta);
};

#endif