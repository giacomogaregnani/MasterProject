#ifndef RAM_HPP
#define RAM_HPP

#include <Eigen/Dense>

using namespace Eigen;

class RAMParam {
private:

    // RAM matrix
    MatrixXd S;

    // Desired acceptance ratio
    double targetAlpha;

    // Counter
    unsigned long count;

    // Size of the parameter space
    unsigned int nParam;

public:

    RAMParam() {};

    void init(double gamma, double alpha, unsigned int paramSize);

    void update(VectorXd& w, double alpha);

    MatrixXd& getS(void);
};

#endif