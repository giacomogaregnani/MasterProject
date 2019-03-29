#include "RandomTimeStep.hpp"

RungeKuttaRandomH::RungeKuttaRandomH(std::default_random_engine* generator,
                                     odeDef ODE,
                                     Butcher tableau,
                                     double meanTimeStep,
                                     double orderP,
                                     double C):
        noiseConstant(C),
        ODE(ODE)
{
    randomGenerator = generator;
    currentTime = 0.0;
    RKMethod = std::make_shared<RungeKutta>(ODE, tableau);
    meanStep = meanTimeStep;
    timeStepOrder = orderP;

    std::uniform_real_distribution<double>::param_type distrParam(meanStep - noiseConstant * std::pow(meanStep, timeStepOrder),
                                                                  meanStep + noiseConstant * std::pow(meanStep, timeStepOrder));
    unifDistribution.param(distrParam);

}

VectorXd RungeKuttaRandomH::oneStep(VectorXd& lastValue, VectorXd& parameter)
{
    double h = unifDistribution(*randomGenerator);
    // double h = lognormalDistribution(*randomGenerator);
    currentTime += h;
    return RKMethod->oneStep(h, lastValue, parameter);
}

void RungeKuttaRandomH::setH(double h)
{
    meanStep = h;
    std::uniform_real_distribution<double>::param_type distrParam(meanStep - noiseConstant * std::pow(meanStep, timeStepOrder),
                                                                  meanStep + noiseConstant * std::pow(meanStep, timeStepOrder));
    unifDistribution.param(distrParam);
}

void RungeKuttaRandomH::setConstant(double C)
{
    noiseConstant = C;
    std::uniform_real_distribution<double>::param_type distrParam(meanStep - noiseConstant * std::pow(meanStep, timeStepOrder),
                                                                  meanStep + noiseConstant * std::pow(meanStep, timeStepOrder));
    unifDistribution.param(distrParam);
}