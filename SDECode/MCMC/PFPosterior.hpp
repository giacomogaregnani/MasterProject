#ifndef PFPOSTERIOR_HPP
#define PFPOSTERIOR_HPP

#include "Posterior.hpp"
#include <EulerMaruyama.hpp>
#include <ParFil.hpp>
#include <Eigen/Dense>
#include <vector>
#include <memory>

// This class computes the posterior for the parameters of an SDE with the Boostrap Particle Filter

using namespace Eigen;

class PFPosterior : public Posterior {
private:
    std::shared_ptr<ParFil> ParticleFilter;
    bool IS;

public:
    PFPosterior() = default;
    virtual ~PFPosterior() = default;
    PFPosterior(std::vector<double>& x, double T, double IC,
                unsigned int sR, double noise, oneDimSde sde,
                double eps, unsigned long M, bool IS = false,
                std::vector<double> timeNoise = {});
    double computePosterior(VectorXd& theta) override;
};

#endif