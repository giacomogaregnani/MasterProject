#ifndef RANDOMCALIBRATION_HPP
#define RANDOMCALIBRATION_HPP

#include "RandomTimeStep.hpp"
#include "AdditiveNoise.hpp"

template<class T> class RandomCalibration {
private:
    double finalTime;

    T* probMethod;

    RungeKutta* detMethod;

    RungeKutta* embMethod;

    std::vector<VectorXd> detSolution;

    std::vector<VectorXd> embErrors;

    void embErrorsEstimate();

    double distribution(double C, unsigned int nMC);

    double h;

    odeDef ODE;

    std::vector<double> samplesConstant;

    std::vector<double> densitiesConstant;

public:
    RandomCalibration()  {};

    RandomCalibration(double finalTime, T* method,
                      RungeKutta* detMethod, RungeKutta* embMethod);

    void calibrate(unsigned int nMC, unsigned int nMCMC, double sigma);

    std::vector<double>& getConstants() {return samplesConstant;};

    std::vector<double>& getDensities() {return densitiesConstant;};
};


#endif