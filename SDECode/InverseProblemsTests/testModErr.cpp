#include <fstream>
#include <iostream>
#include <iomanip>
#include "generateObservations.hpp"
#include <MCMC.hpp>

// Remark: epsilon = p(0)

double gradV0(double x)
{
    return x;
}

double V1(double x)
{
    return std::cos(x);
}

double gradV1(double x)
{
    return -std::sin(x);
}

double multiDrift(double x, VectorXd& p)
{
    return -1.0 * (std::exp(p(1)) * gradV0(x) + 1 / p(0) * gradV1(x / p(0)));
}

double homoDrift(double x, VectorXd& p)
{
    return -1.0 * std::exp(p(1)) * gradV0(x);
}

double diffusion(double x, VectorXd &p)
{
    return std::sqrt(2.0 * std::exp(p(2)));
}

int main(int argc, char* argv[])
{
    oneDimSde sde;
    sde.drift = &multiDrift;
    sde.diffusion = &diffusion;
    oneDimSde sdeHomo;
    sdeHomo.drift = &homoDrift;
    sdeHomo.diffusion = &diffusion;

    double T = 1;
    unsigned int N = 10000;

    VectorXd param(3);
    param(0) = 0.01;  // Epsilon
    param(1) = std::log(1.0);  // True multiscale alpha
    param(2) = std::log(0.5);  // True multiscale betainv

    std::ofstream output(DATA_PATH + std::string("testModErr.txt"), std::ofstream::out | std::ofstream::trunc);

    // Modeling error estimation
    double IC = 0.0;
    VectorXd priorMean(param.size());
    priorMean << param(0), 0.0, 0.0;
    VectorXd priorStdDev(param.size());
    priorStdDev << 0.0, 1.0, 1.0;

    // Generate and perturb observations
    double noise = 1e-2;
    std::default_random_engine noiseSeed{1};
    std::normal_distribution<double> noiseDistribution(0.0, noise);
    std::cout << "Generating observations..." << std::endl;
    auto x = generateObservations1D(sde, IC, param, T, N, 0);
    std::cout << "Generated observations" << std::endl;
    for (unsigned long i = 1; i < x.size(); i++) {
        x[i] += noiseDistribution(noiseSeed);
    }

    ModErr modErr(sdeHomo, sde, &V1, IC, priorMean, priorStdDev, T, N, x, 0);
    unsigned int nMC = 400;
    unsigned int nParam = 10;

    modErr.compute(nParam, nMC);
    std::vector<double> means, stddevs;
    modErr.getStats(means, stddevs);

    for (unsigned int i = 0; i < N+1; i++) {
        output << means[i] << "\t";
    }
    output << std::endl;
    for (unsigned int i = 0; i < N+1; i++) {
        output << stddevs[i] << "\t";
    }
    output << std::endl;

    modErr.computePF(nParam, nMC);
    modErr.getStats(means, stddevs);

    for (unsigned int i = 0; i < N+1; i++) {
        output << means[i] << "\t";
    }
    output << std::endl;

    for (unsigned int i = 0; i < N+1; i++) {
        output << stddevs[i] << "\t";
    }
    output << std::endl;

    output.close();
    return 0;
}