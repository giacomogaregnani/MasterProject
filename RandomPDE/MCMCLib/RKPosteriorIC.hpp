#ifndef RKPOSTERIORIC_HPP
#define RKPOSTERIORIC_HPP

#include <Eigen/Dense>
#include <RungeKuttaSolver.hpp>
#include "Proposals.hpp"
#include "Posterior.hpp"
#include "KarhunenLoeve.hpp"

using namespace Eigen;

class RKPosteriorIC : public Posterior {
private:

    odeDef ODE;

    RungeKutta RKSolver;

    double h;

    std::vector<VectorXd> observations;

    std::vector<double> tObs;

    double noise;

    MatrixXd noiseMatrix;

    LLT<MatrixXd> noiseMatrixLLT;

    bool isIdentity;

public:

    RKPosteriorIC() {};

    RKPosteriorIC(double h,
                  std::vector<VectorXd>& observations,
                  std::vector<double>& tObs,
                  double noise,
                  odeDef ODE,
                  Butcher tableau);

    RKPosteriorIC(double h,
                  std::vector<VectorXd>& observations,
                  std::vector<double>& tObs,
                  MatrixXd noise,
                  odeDef ODE,
                  Butcher tableau);

    double computePosterior(VectorXd& theta);

    VectorXd coeffToField(VectorXd& theta) {};
};


#endif