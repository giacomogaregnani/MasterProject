#ifndef KARHUNENLOEVE_HPP
#define KARHUNENLOEVE_HPP

#include <Eigen/Dense>
#include <vector>

#ifndef PI
#define PI 3.14159265359
#endif

using namespace Eigen;

enum covKernelFcts {
    INVLAPDIR,
    SQDEXP,
    EXP
};

class KarhunenLoeve {
private:
    VectorXd mean;

    MatrixXd covariance;

    void buildCovariance(covKernelFcts cov);

    std::vector<double> param;

    std::vector<VectorXd> covEVec;

    std::vector<double> covSqrtEVal;

    int nElem;

public:

    KarhunenLoeve() = default;

    KarhunenLoeve(VectorXd& fieldMean, covKernelFcts cov, int N, std::vector<double> param = {});

    VectorXd KL(VectorXd& coeffiecients);

    std::vector<VectorXd> getEVecs(void);
};

#endif
