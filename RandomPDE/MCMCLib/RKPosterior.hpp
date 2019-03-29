#ifndef RKPOSTERIOR_HPP
#define RKPOSTERIOR_HPP

#include <Eigen/Dense>
#include <RungeKuttaSolver.hpp>
#include "Proposals.hpp"
#include "Posterior.hpp"

using namespace Eigen;

class RKPosterior : public Posterior {
private:
    RungeKutta RKSolver;

    double h;

    double T;

    VectorXd IC;

    std::vector<VectorXd> observations;

    std::vector<double> tObs;

    double noise;

public:
    RKPosterior() {};

    RKPosterior(double h,
                double T,
                VectorXd& initialCondition,
                std::vector<VectorXd>& observations,
                std::vector<double>& tObs,
                double noise,
                odeDef ODE,
                Butcher tableau);

    double computePosterior(VectorXd& theta);

    VectorXd coeffToField(VectorXd& theta) {return VectorXd::Zero(0);}
};

#endif
