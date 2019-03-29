#ifndef ADDNOISE_HPP
#define ADDNOISE_HPP

#include <RungeKuttaSolver.hpp>
#include <memory>
#include <random>

class RungeKuttaAddNoise {

private:
    std::shared_ptr<RungeKutta> RKMethod;

    std::default_random_engine* randomGenerator;

    std::normal_distribution<double> normalDistribution;

    double timeStep;

    double noiseOrder;

    double noiseConstant;

    odeDef ODE;

public:
    VectorXd currentNoise;

    RungeKuttaAddNoise() {};

    RungeKuttaAddNoise(std::default_random_engine* generator,
                       odeDef ODE,
                       Butcher tableau,
                       double meanTimeStep,
                       double p,
                       double C = 1);

    VectorXd oneStep(VectorXd& lastValue, VectorXd& parameter);

    void setH(double h);

    void setConstant(double C);

    double getH(void) {return timeStep;};

    odeDef getODE(void) {return ODE;};
};

#endif