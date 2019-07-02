#ifndef PFPOSTERIORHOM_HPP
#define PFPOSTERIORHOM_HPP

#include "Posterior.hpp"
#include <EulerMaruyama.hpp>
#include <ParFil.hpp>
#include <Eigen/Dense>
#include <vector>
#include <memory>

// This class computes the posterior for the parameters of an SDE with the Particle Filter

using namespace Eigen;

class PFPosteriorHom : public Posterior {
private:
    std::shared_ptr<ParFil> ParticleFilter;
    bool IS;
    double (*V1) (double);
    double eps;
    std::vector<std::vector<double>>* errors;

public:
    PFPosteriorHom() = default;
    virtual ~PFPosteriorHom() = default;
    VectorXd computeHomogeneous(VectorXd param, double L, double (*V1) (double));
    PFPosteriorHom(std::vector<double>& x, double T, double IC,
                   unsigned int sR, double noise, oneDimSde sde,
                   double (*V1) (double),
                   double eps, unsigned long M, bool IS = false,
                   std::vector<double> timeNoise = {},
                   std::vector<std::vector<double>>* errors = nullptr);
    double computePosterior(VectorXd& theta) override;
};

#endif