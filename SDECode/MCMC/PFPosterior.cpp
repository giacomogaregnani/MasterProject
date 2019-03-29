#include <iostream>
#include "PFPosterior.hpp"

PFPosterior::PFPosterior(std::vector<double>& x, double T, double IC,
                         unsigned int sR, double noise, oneDimSde sde,
                         double eps, unsigned long M)
{
    ParticleFilter = std::make_shared<ParFil>(x, T, IC, sR, noise, sde, eps, M);
}

double PFPosterior::computePosterior(VectorXd theta)
{
    VectorXd thetaWithoutEps(theta.size()-1);
    thetaWithoutEps = theta.tail(theta.size()-1);
    double prior = -0.5 * thetaWithoutEps.dot(thetaWithoutEps);

    ParticleFilter->compute(theta);
    double likelihood = ParticleFilter->getLikelihood();

    return prior + likelihood;
}