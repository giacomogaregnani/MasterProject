#ifndef RKPROBPOSTERIOR_HPP
#define RKPROBPOSTERIOR_HPP

#include <Eigen/Dense>
#include <RandomTimeStep.hpp>
#include "Proposals.hpp"
#include "Posterior.hpp"

using namespace Eigen;

class RKProbPosterior : public Posterior {
private:
    std::vector<RungeKuttaAddNoise> RKSolvers;

    double h;

    double T;

    VectorXd IC;

    std::vector<VectorXd> observations;

    std::vector<double> tObs;

    double noise;

    int nMC;

    double p;

    std::vector<std::default_random_engine> generators;

public:
    RKProbPosterior() {};

    RKProbPosterior(double h,
                    double T,
                    VectorXd& initialCondition,
                    std::vector<VectorXd>& observations,
                    std::vector<double>& tObs,
                    double noise,
                    odeDef ODE,
                    Butcher tableau,
                    int nMC, double p);

    double computePosterior(VectorXd& theta);

    VectorXd coeffToField(VectorXd& theta) {return VectorXd::Zero(0);}
};

#endif
