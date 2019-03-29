#include "AdditiveNoise.hpp"

RungeKuttaAddNoise::RungeKuttaAddNoise(std::default_random_engine *generator,
                                       odeDef ODE,
                                       Butcher tableau, double h,
                                       double p,
                                       double C):
    noiseConstant(C),
    ODE(ODE)
{
    randomGenerator = generator;
    RKMethod = std::make_shared<RungeKutta>(ODE, tableau);
    timeStep = h;
    noiseOrder = p;
    std::normal_distribution<double>::param_type distrParam(0.0, noiseConstant * std::pow(timeStep, noiseOrder + 0.5));
    normalDistribution.param(distrParam);
    currentNoise.resize(ODE.size);
}

VectorXd RungeKuttaAddNoise::oneStep(VectorXd& lastValue, VectorXd& parameter)
{
    for (int i = 0; i < ODE.size; i++) {
        currentNoise(i) = normalDistribution(*randomGenerator);
    }

    return RKMethod->oneStep(timeStep, lastValue, parameter) + currentNoise;
}

void RungeKuttaAddNoise::setH(double h) {
    timeStep = h;
    std::normal_distribution<double>::param_type distrParam(0.0, noiseConstant * std::pow(timeStep, noiseOrder + 0.5));
    normalDistribution.param(distrParam);
    std::cout << std::endl;
}

void RungeKuttaAddNoise::setConstant(double C)
{
    noiseConstant = C;
    std::normal_distribution<double>::param_type distrParam(0.0, noiseConstant * std::pow(timeStep, noiseOrder + 0.5));
    normalDistribution.param(distrParam);
}