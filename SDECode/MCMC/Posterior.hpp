#ifndef POSTERIOR_HPP
#define POSTERIOR_HPP

#include <Eigen/Dense>

using namespace Eigen;

class Posterior {
public:
    virtual ~Posterior() = default;
    virtual double computePosterior(VectorXd& theta) = 0;
};

#endif
