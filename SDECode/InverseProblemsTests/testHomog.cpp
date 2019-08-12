#include "generateObservations.hpp"
#include "computeHomogeneous.hpp"
#include "../matplotlib-cpp-master/matplotlibcpp.h"

namespace plt = matplotlibcpp;

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
    oneDimSde sde{&multiDrift, &diffusion};
    oneDimSde sdeHomo{&homoDrift, &diffusion};

    double T = 1;
    unsigned int N = 40000;

    VectorXd tmpParam(3);
    tmpParam(0) = 0.025;  // Epsilon
    tmpParam(1) = 1.0;  // True multiscale alpha
    tmpParam(2) = 0.5;  // True multiscale betainv

    // Compute coefficients of the homogenised equation
    std::vector<double> homCoeffs = computeHomCoeffs(tmpParam, (2.0*M_PI), &V1);
    VectorXd homoParam(3);
    homoParam(0) = 0.1;
    homoParam(1) = std::log(homCoeffs[0]);
    homoParam(2) = std::log(homCoeffs[1]);
    VectorXd param = tmpParam;
    param(1) = std::log(param(1));
    param(2) = std::log(param(2));

    double IC = 0.0;
    std::random_device dev;
    unsigned int obsSeed;

    unsigned long M = 5000;
    std::vector<double> xFinal(M), xHomFinal(M);

    std::vector<double> timeVec(N+1);
    for (unsigned int k = 0; k < N+1; k++) {
        timeVec[k] = k*T/N;
    }

    #pragma omp parallel for num_threads(5)
    for (unsigned int k = 0; k < M; k++) {
        obsSeed = dev();
        auto x = generateObservations1D(sde, IC, param, T, N, obsSeed);
        xFinal[k] = x.back();

        // Compute a homogeneous path with the same Brownian for comparison
        auto xHom = generateObservations1D(sdeHomo, IC, homoParam, T, N, obsSeed);
        xHomFinal[k] = xHom.back();

        if (k < 1) {
            plt::plot(timeVec, xHom, "b");
            plt::plot(timeVec, x, "r");
        }
    }
    plt::show();

    plt::hist(xHomFinal, 20, "b", 0.3);
    plt::hist(xFinal, 20, "r", 0.3);
    plt::show();

    return 0;
}