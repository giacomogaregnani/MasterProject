#include "EulerMaruyama.hpp"
#include <iostream>


EM::EM(multiDimSde& sde, VectorXd& param, std::default_random_engine& seed):
        sde(sde),
        param(param),
        seed(seed)
{
    dBM = std::normal_distribution<double>(0, 1);
}

VectorXd EM::oneStep(double h, VectorXd Xn)
{
    VectorXd dW(sde.nBM);
    for (unsigned int i = 0; i < sde.nBM; i++)
        dW(i) = dBM(seed);
    return Xn + h * sde.drift(Xn, param) + sde.diffusion(Xn, param) * std::sqrt(h) * dW;
}

void EM::modifyParam(VectorXd& newParam)
{
    param = newParam;
}

// ONE DIMENSIONAL SOLVER

EM1D::EM1D(oneDimSde& sde, VectorXd& param, std::default_random_engine& seed):
        sde(sde),
        param(param),
        seed(seed)
{
    dBM = std::normal_distribution<double>(0, 1);
}

EM1D::EM1D(oneDimSde& sde, std::default_random_engine& seedIn):
        sde(sde),
        seed(seed)
{
    param = VectorXd::Zero(1);
    dBM = std::normal_distribution<double>(0, 1);
}

double EM1D::oneStep(double h, double Xn)
{
    return Xn + h * sde.drift(Xn, param) + sde.diffusion(Xn, param) * std::sqrt(h) * dBM(seed);
}

double EM1D::oneStepGivenNoise(double h, double Xn, double noise)
{
    return Xn + h * sde.drift(Xn, param) + sde.diffusion(Xn, param) * noise;
}


void EM1D::modifyParam(VectorXd& newParam)
{
    param = newParam;
}