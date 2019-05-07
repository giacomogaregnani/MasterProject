#ifndef RKPROBPOSTERIORIC_HPP
#define RKPROBPOSTERIORIC_HPP

#include <Eigen/Dense>
#include <RandomTimeStep.hpp>
#include "Posterior.hpp"

using namespace Eigen;

class RKProbPosteriorIC : public Posterior {
private:
    odeDef ODE;
    std::vector<RungeKuttaRandomH> RKSolvers;
    double h;
    std::vector<VectorXd> observations;
    std::vector<double> tObs;
    double noise;
    int nMC;
    double p;
    std::vector<std::default_random_engine> generators;

public:
    RKProbPosteriorIC() {};
    RKProbPosteriorIC(double h,
                      std::vector<VectorXd>& observations,
                      std::vector<double>& tObs,
                      double noise,
                      odeDef ODE,
                      Butcher tableau,
                      int nMC, double p);
    double computePosterior(VectorXd& theta);
    VectorXd coeffToField(VectorXd& theta) {};
};


#endif