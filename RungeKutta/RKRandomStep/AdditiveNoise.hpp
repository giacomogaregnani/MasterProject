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

    int ODESize;

public:
    VectorXd currentNoise;

    RungeKuttaAddNoise(std::default_random_engine* generator,
                       odeDef ODE, std::vector<double> paramVec,
                       Butcher tableau,
                       double meanTimeStep,
                       double orderP);

    VectorXd oneStep(VectorXd lastValue);
};