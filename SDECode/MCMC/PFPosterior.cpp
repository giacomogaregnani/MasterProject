#include <iostream>
#include "PFPosterior.hpp"

PFPosterior::PFPosterior(std::vector<double>& x, double T, double IC,
                         double noise, oneDimSde sde,
                         double eps, unsigned long M, bool IS):
                         IS(IS)
{
    std::shared_ptr<ForwardPF> forwardSampler;
    auto N = x.size() - 1;
    std::random_device dev;
    std::default_random_engine seed(dev());
    forwardSampler = std::make_shared<ForwardPF>(sde, T/N, seed);
    ParticleFilter = std::make_shared<ParFil>(x, T, IC, noise, eps, M, forwardSampler);
}

double PFPosterior::computePosterior(VectorXd& theta)
{
    VectorXd thetaWithoutEps(theta.size()-1);
    thetaWithoutEps = theta.tail(theta.size()-1);
    double prior = -0.5 * thetaWithoutEps.dot(thetaWithoutEps);

    if (IS) {
        ParticleFilter->computeDiffBridge(theta);
    } else {
        ParticleFilter->compute(theta);
    }
    double likelihood = ParticleFilter->getLikelihood();

    return prior + likelihood;
}