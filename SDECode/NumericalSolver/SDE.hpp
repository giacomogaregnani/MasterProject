#include <Eigen/Dense>

using namespace Eigen;

struct multiDimSde {
    unsigned int size;
    unsigned int nBM;
    VectorXd (*drift) (VectorXd, VectorXd&);
    MatrixXd (*diffusion) (VectorXd, VectorXd&);
};

struct oneDimSde {
    double (*drift) (double, VectorXd&);
    double (*diffusion) (double, VectorXd&);
};