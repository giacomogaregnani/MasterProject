#include <iostream>
#include "PFPosteriorHom.hpp"

PFPosteriorHom::PFPosteriorHom(std::vector<double>& x, double T, double IC,
                               unsigned int sR, double noise, oneDimSde sde,
                               double (*V1) (double),
                               double eps, unsigned long M, bool IS,
                               std::vector<double> timeNoise):
        IS(IS),
        V1(V1),
        eps(eps)
{
    ParticleFilter = std::make_shared<ParFil>(x, T, IC, sR, noise, sde, eps, M, timeNoise);
}

VectorXd PFPosteriorHom::computeHomogeneous(VectorXd param, double L, double (*V1) (double))
{
    param(1) = std::exp(param(1));
    param(2) = std::exp(param(2));

    // Compute the coefficients Z and \hat{Z}

    unsigned int N = 10000;
    VectorXd discr;
    discr.setLinSpaced(N+1, 0.0, L);
    double h = L / N;

    std::vector<double> Zs = {0.0, 0.0};

    for (int i = 0; i < N; i++) {
        Zs[0] += h * (std::exp( V1(discr[i]) / param(2)) + std::exp( V1(discr[i+1]) / param(2))) / 2;
        Zs[1] += h * (std::exp(-V1(discr[i]) / param(2)) + std::exp(-V1(discr[i+1]) / param(2))) / 2;
    }

    VectorXd homParam(2);

    homParam(0) = param(1) * L * L / (Zs[0] * Zs[1]);
    homParam(1) = param(2) * L * L / (Zs[0] * Zs[1]);
    homParam(0) = std::log(homParam(0));
    homParam(1) = std::log(homParam(1));

    return homParam;
}

double PFPosteriorHom::computePosterior(VectorXd& theta)
{
    VectorXd thetaWithoutEps(theta.size()-1);
    thetaWithoutEps = theta.tail(theta.size()-1);
    double prior = -0.5 * thetaWithoutEps.dot(thetaWithoutEps);
    theta(0) = eps;

    VectorXd thetaHom(3);
    VectorXd tmp = computeHomogeneous(theta, (2.0*M_PI), V1);

    thetaHom(0) = theta(0);
    thetaHom(1) = tmp(0);
    thetaHom(2) = tmp(1);

    if (IS) {
        ParticleFilter->computeDiffBridge(thetaHom);
    } else {
        ParticleFilter->compute(thetaHom);
    }
    double likelihood = ParticleFilter->getLikelihood();

    return prior + likelihood;
}
