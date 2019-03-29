#ifndef RANDOMSTEP_HPP
#define RANDOMSTEP_HPP

#include "AdditiveNoise.hpp"

class RungeKuttaRandomH {

private:
    std::shared_ptr<RungeKutta> RKMethod;

    std::default_random_engine* randomGenerator;

    std::uniform_real_distribution<double> unifDistribution;

    double currentTime;

    double meanStep;

    double timeStepOrder;

    double noiseConstant;

    odeDef ODE;

public:
    RungeKuttaRandomH() {};

    RungeKuttaRandomH(std::default_random_engine *generator,
                      odeDef ODE,
                      Butcher tableau,
                      double meanTimeStep,
                      double orderP,
                      double C = 1);

    VectorXd oneStep(VectorXd& lastValue, VectorXd& parameter);

    void setH(double h);

    void setConstant(double C);

    double getH(void) {return meanStep;};

    odeDef getODE(void) {return ODE;}
};

#endif