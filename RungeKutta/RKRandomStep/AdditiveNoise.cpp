#include "AdditiveNoise.hpp"

RungeKuttaAddNoise::RungeKuttaAddNoise(std::default_random_engine *generator,
                                       odeDef ODE, std::vector<double> paramVec,
                                       Butcher tableau, double h,
                                       double orderP)
{
    randomGenerator = generator;
    RKMethod = std::make_shared<RungeKutta>(ODE, paramVec, tableau);
    timeStep = h;
    noiseOrder = orderP;
    std::normal_distribution<double>::param_type distrParam(0.0,
                                                            std::pow(timeStep, noiseOrder + 0.5));
    normalDistribution.param(distrParam);
    ODESize = ODE.size;
    currentNoise.resize(ODESize);
}

VectorXd RungeKuttaAddNoise::oneStep(VectorXd lastValue)
{
    for (int i = 0; i < ODESize; i++) {
        currentNoise(i) = normalDistribution(*randomGenerator);
    }

    return RKMethod->oneStep(lastValue, timeStep) + currentNoise;
}
