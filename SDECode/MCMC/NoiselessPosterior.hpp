#ifndef NOISELESSPOSTERIOR_HPP
#define NOISELESSPOSTERIOR_HPP

#include "Posterior.hpp"
#include <EulerMaruyama.hpp>
#include <ParFilLib.hpp>
#include <Eigen/Dense>
#include <vector>
#include <memory>

using namespace Eigen;

class NLPosterior : public Posterior {
private:
    std::vector<double> obs;
    double h;
    double (*gradV0) (double);
    VectorXd paramSde;
    unsigned long N;

public:
    NLPosterior() = default;
    virtual ~NLPosterior() = default;
    NLPosterior(std::vector<double>& x, double h,
                double (*gradV0) (double), VectorXd& param);
    double computePosterior(VectorXd& theta) override;
};


#endif
