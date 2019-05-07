#ifndef RKADDPOSTERIORIC_HPP
#define RKADDPOSTERIORIC_HPP

#include <Eigen/Dense>
#include <RandomTimeStep.hpp>
#include "Posterior.hpp"

using namespace Eigen;

class RKAddPosteriorIC : public Posterior {
private:
    odeDef ODE;
    std::vector<RungeKuttaAddNoise> RKSolvers;
    double h;
    std::vector<VectorXd> observations;
    std::vector<double> tObs;
    double noise;
    int nMC;
    double p;
    std::vector<std::default_random_engine> generators;

public:
    RKAddPosteriorIC() {};
    RKAddPosteriorIC(double h,
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