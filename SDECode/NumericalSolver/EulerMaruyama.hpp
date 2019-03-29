#ifndef EULERMARUYAMA_HPP
#define EULERMARUYAMA_HPP

#include "SDE.hpp"
#include <random>
#include <iostream>

// Euler-Maruyama in the multi-d case

class EM {
private:
    multiDimSde sde;
    VectorXd param;
    std::default_random_engine seed;
    std::normal_distribution<double> dBM;

public:
    EM() = default;
    EM(multiDimSde& sde, VectorXd& param, std::default_random_engine& seed);
    VectorXd oneStep(double h, VectorXd Xn);
    void modifyParam(VectorXd& newParam);
};

// Euler-Maruyama in the 1d case

class EM1D {
    oneDimSde sde;
    VectorXd param;
    std::default_random_engine seed;
    std::normal_distribution<double> dBM;

public:
    EM1D() = default;
    ~EM1D() = default;
    EM1D(oneDimSde& sde, std::default_random_engine& seed);
    EM1D(oneDimSde& sde, VectorXd& param, std::default_random_engine& seed);
    double oneStep(double h, double Xn);
    double oneStepGivenNoise(double h, double Xn, double noise);
    void modifyParam(VectorXd& newParam);
};

#endif