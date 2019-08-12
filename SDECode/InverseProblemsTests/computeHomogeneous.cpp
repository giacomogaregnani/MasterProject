#include "computeHomogeneous.hpp"


std::vector<double> computeHomCoeffs(VectorXd& param, double L, double (*V1)(double))
{
    std::vector<double> homParam(2);
    auto Zs = computeZs(V1, param(2), L);

    homParam[0] = param(1) * L * L / (Zs[0] * Zs[1]);
    homParam[1] = param(2) * L * L / (Zs[0] * Zs[1]);

    return homParam;
}

std::vector<double> computeZs(double (*gradV1)(double), double sigma, double L)
{
    unsigned int N = 10000;
    VectorXd discr;
    discr.setLinSpaced(N+1, 0.0, L);
    double h = L / N;

    std::vector<double> Zs = {0.0, 0.0};

    for (int i = 0; i < N; i++) {
        Zs[0] += h * (std::exp( gradV1(discr[i]) / sigma) + std::exp( gradV1(discr[i+1]) / sigma)) / 2;
        Zs[1] += h * (std::exp(-gradV1(discr[i]) / sigma) + std::exp(-gradV1(discr[i+1]) / sigma)) / 2;
    }

    return Zs;
}
