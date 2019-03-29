#include "KarhunenLoeve.hpp"
#include <iostream>

KarhunenLoeve::KarhunenLoeve(VectorXd &fieldMean, covKernelFcts cov,
                             int N, std::vector<double> param):
    param(param),
    nElem(N - 1)
{
    mean = fieldMean;
    covariance.resize(nElem + 1, nElem + 1);
    covEVec = std::vector<VectorXd>(nElem + 1, VectorXd(nElem + 1));
    covSqrtEVal = std::vector<double>(nElem + 1, 0.0);
    buildCovariance(cov);
}

void KarhunenLoeve::buildCovariance(covKernelFcts cov) {

    double h = 1.0 / nElem;
    MatrixXd EVec;
    VectorXd EVal;
    EigenSolver<MatrixXd> es;

    switch (cov) {
        case INVLAPDIR:
            for (int j = 0; j < nElem + 1; j++) {
                covSqrtEVal[j] = 1.0 / (PI * (j + 1));

                covEVec[j](0) = 0;
                covEVec[j](nElem) = 0;
                for (int i = 1; i < nElem; i++) {
                    covEVec[j](i) = sin(PI * (j + 1) * (h * i));
                }
                covEVec[j] *= sqrt(2.0);
            }
            break;

        case SQDEXP:
            for (int i = 0; i < nElem + 1; i++) {
                for (int j = 0; j < nElem + 1; j++) {
                    covariance(i, j) = param[0] * std::exp(-(h * i - h * j) * (h * i - h * j) / (2 * param[1] * param[1]));
                }
            }

            es = EigenSolver<MatrixXd>(covariance);
            EVec = es.eigenvectors().real();
            EVal = es.eigenvalues().real();

            for (int i = 0; i < nElem + 1; i++) {
                covEVec[i] = EVec.col(i);
                covSqrtEVal[i] = EVal(i);
            }
            break;

        case EXP:
            for (int i = 0; i < nElem + 1; i++) {
                for (int j = 0; j < nElem + 1; j++) {
                    covariance(i, j) = param[0] * std::exp(-std::abs(h * i - h * j) / (2 * param[1] * param[1]));
                }
            }

            es = EigenSolver<MatrixXd>(covariance);
            EVec = es.eigenvectors().real();
            EVal = es.eigenvalues().real();

            for (int i = 0; i < nElem + 1; i++) {
                covEVec[i] = EVec.col(i);
                covSqrtEVal[i] = EVal(i);
            }
            break;
    }
}

VectorXd KarhunenLoeve::KL(VectorXd& coeffs)
{
    VectorXd sample = mean;

    for (int i = 0; i < coeffs.size(); i++) {
        sample += covEVec[i] * covSqrtEVal[i] * coeffs(i);
    }

    return sample;
}

std::vector<VectorXd> KarhunenLoeve::getEVecs()
{
    return covEVec;
}
