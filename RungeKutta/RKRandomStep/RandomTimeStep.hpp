#include "AdditiveNoise.hpp"

class RungeKuttaRandomH {

private:
    std::shared_ptr<RungeKutta> RKMethod;

    std::default_random_engine* randomGenerator;

    std::uniform_real_distribution<double> unifDistribution;

    double currentTime;

    double meanStep;

    double timeStepOrder;

public:
    RungeKuttaRandomH(std::default_random_engine* generator,
                      odeDef ODE, std::vector<double> paramVec,
                      Butcher tableau,
                      double meanTimeStep,
                      double orderP);

    VectorXd oneStep(VectorXd lastValue);

    double getCurrentTime(void);

    void resetTimeToZero();
};