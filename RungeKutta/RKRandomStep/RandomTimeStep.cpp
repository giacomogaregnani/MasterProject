#include "RandomTimeStep.hpp"

RungeKuttaRandomH::RungeKuttaRandomH(std::default_random_engine* generator,
                                     odeDef ODE,
                                     std::vector<double> paramVec,
                                     Butcher tableau,
                                     double meanTimeStep,
                                     double orderP)
{
    randomGenerator = generator;
    currentTime = 0.0;
    RKMethod = std::make_shared<RungeKutta>(ODE, paramVec, tableau);
    meanStep = meanTimeStep;
    timeStepOrder = orderP;

    std::uniform_real_distribution<double>::param_type distrParam(meanStep - std::pow(meanStep, timeStepOrder),
                                                                  meanStep + std::pow(meanStep, timeStepOrder));

    unifDistribution.param(distrParam);
}

VectorXd RungeKuttaRandomH::oneStep(VectorXd lastValue)
{
    double h = unifDistribution(*randomGenerator);
    currentTime += h;
    return RKMethod->oneStep(lastValue, h);
}

double RungeKuttaRandomH::getCurrentTime(void)
{
    return currentTime;
}

void RungeKuttaRandomH::resetTimeToZero()
{
    currentTime = 0.0;
}