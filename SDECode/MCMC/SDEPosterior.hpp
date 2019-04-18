#ifndef SDEPOSTERIOR_HPP
#define SDEPOSTERIOR_HPP

#include "Posterior.hpp"
#include <EulerMaruyama.hpp>
#include <Eigen/Dense>
#include <vector>

using namespace Eigen;

class SDEPosterior : public Posterior {
private:
    EM1D Solver;
    double T;
    double IC;
    unsigned int samplingRatio;
    std::vector<double> obs;
    double noise;
    unsigned long nMC;
    double eps;

public:
    SDEPosterior() = default;
    virtual ~SDEPosterior() = default;
    SDEPosterior(std::vector<double>& x, double T, double IC,
                 unsigned int sR, double noise, oneDimSde sde,
                 double eps, unsigned long M);
    double computePosterior(VectorXd& theta) override;
};

#endif