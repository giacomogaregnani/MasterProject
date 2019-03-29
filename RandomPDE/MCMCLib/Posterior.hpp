#ifndef POSTERIOR_HPP
#define POSTERIOR_HPP

#include <Eigen/Dense>

using namespace Eigen;

class Posterior {
public:
    virtual double computePosterior(VectorXd& theta) = 0;

    virtual VectorXd coeffToField(VectorXd& theta) = 0;
};

#endif
