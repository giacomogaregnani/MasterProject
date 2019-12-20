#include <iostream>
#include "generateObservations.hpp"
#include "Filter.hpp"
#include "../matplotlib-cpp-master/matplotlibcpp.h"
#include "computeEstimators.hpp"
#include "computeHomogeneous.hpp"

// Remark: epsilon = p(0)

namespace plt = matplotlibcpp;

double gradV0(double x)
{
    return x;
}

double laplV0(double x)
{
    return 1.0;
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
    return -1.0 * (p(1) * gradV0(x) + 1 / p(0) * gradV1(x / p(0)));
}

double diffusion(double x, VectorXd &p)
{
    return std::sqrt(2.0 * p(2));
}

int main(int argc, char* argv[])
{
    oneDimSde sde{&multiDrift, &diffusion};
    double IC = 0.0;

    VectorXd param(3);
    param(0) = 0.05;
    param(1) = 1.0; // True multiscale alpha
    param(2) = 0.7; // True multiscale sigma

    std::vector<double> avg;

    // Refer to the Caltech notes
    double beta = 2.5;
    double zeta = 0.5;
    double gamma = -std::log(60.) / std::log(param(0));

    std::random_device dev;
    std::default_random_engine proposalSeed(dev());
    std::default_random_engine acceptanceSeed(dev());

    double eps = param(0);
    std::vector<double> homCoeffs = computeHomCoeffs(param, (2.0*M_PI), &V1);
    std::cout << homCoeffs[0] << " " << homCoeffs[1] << std::endl;

    auto h = std::pow(eps, beta);
    auto N = static_cast<unsigned long>(std::round(std::pow(eps, -beta-gamma)));
    auto T = h * N;

    unsigned int nExp = 1000;

    double betaFilter = 5;
    double deltaFilter = std::pow(eps, zeta);

    std::vector<double> sigmaVec(nExp);
    std::vector<double> aVec(nExp);

    std::vector<double> timeVec(N+1);
    for (unsigned long i = 0; i < N+1; i++) {
        timeVec[i] = h*i;
    }

    #pragma omp parallel for num_threads(10)
    for (unsigned int i = 0; i < nExp; i++) {
        Filter filter(deltaFilter, betaFilter, h, N+1);
        double sigmaCoeff = std::pow(2.0, 1.0/betaFilter) / (h * filter.getFilterZero());
        auto trajectory = generateObservations1D(sde, IC, param, T, N, dev());
        auto filtered = filter.compute(trajectory);
        sigmaVec[i] = sigmaCoeff * estimateSigma(filtered, h);
        aVec[i] = estimateATilde(filtered, h, &gradV0, &laplV0, sigmaVec[i]);
        printf("%d\n", i);
    }

    for (auto const & it : sigmaVec)
        std::cout << it << std::endl;

    // Plots
    std::vector<double> dummyA = {homCoeffs[0], homCoeffs[0]};
    std::vector<double> dummyAVal = {0, 100};

    std::vector<double> dummyS = {homCoeffs[1], homCoeffs[1]};
    std::vector<double> dummySVal = {0, 100};

    plt::hist(aVec, 30, "b", 0.5);
    plt::plot(dummyA, dummyAVal);
    plt::show();

    plt::hist(sigmaVec, 30, "b", 0.5);
    plt::plot(dummyS, dummySVal);
    plt::show();

    return 0;
}