#include "MCMC.hpp"
#include <Eigen/Cholesky>

void RAMParam::init(double gamma, double alpha, unsigned int paramSize)
{
    S = MatrixXd::Identity(paramSize, paramSize) * gamma;
    targetAlpha = alpha;
    count = 1;
    nParam = paramSize;
}

void RAMParam::update(VectorXd& w, double alpha)
{
    MatrixXd WWT = w * w.transpose();
    double WTW = w.dot(w);
    double gammaI = std::min(1.0, 2.0 * pow(static_cast<double>(count++), -0.75));
    double diffAlpha = alpha - targetAlpha;
    double coeff = gammaI * diffAlpha / WTW;
    MatrixXd C = MatrixXd::Identity(nParam, nParam) + WWT * coeff;
    LLT<MatrixXd> chol(S * C * S.transpose());
    S = chol.matrixL();
}

MatrixXd& RAMParam::getS(void)
{
    return S;
}