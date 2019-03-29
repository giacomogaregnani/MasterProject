#ifndef RKPOSTERIORFIXED_HPP
#define RKPOSTERIORFIXED_HPP

#include <Eigen/Dense>
#include <RungeKuttaSolver.hpp>
#include "Proposals.hpp"
#include "Posterior.hpp"

using namespace Eigen;

class RKProbPosteriorFixed : public Posterior {
private:
    RungeKutta RKSolver;

    double h;

    double T;

    VectorXd IC;

    std::vector<VectorXd> observations;

    std::vector<double> tObs;

    double noise;

    std::vector<double> timeSteps;

public:
    RKProbPosteriorFixed() {};

    RKProbPosteriorFixed(double h,
                double T,
                VectorXd& initialCondition,
                std::vector<VectorXd>& observations,
                std::vector<double>& tObs,
                double noise,
                odeDef ODE,
                Butcher tableau,
                double p, int init);

    double computePosterior(VectorXd& theta);

    VectorXd coeffToField(VectorXd& theta) {return VectorXd::Zero(0);}
};

#endif
