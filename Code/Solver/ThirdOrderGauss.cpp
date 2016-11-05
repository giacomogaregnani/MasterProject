#include "Solver.hpp"
#include <unsupported/Eigen/MatrixFunctions>
#include <Eigen/LU>
#include <cmath>

#ifndef PI
#define PI 3.14159265358979323846
#endif

using namespace Eigen;

// Reference to the notation and the formulas in
// "Posterior inference on parameters of Stochastic Differential Equations
// via non-linear Gaussian filtering and adaptive MCMC"
// Simo Särkkä et al., Stat. Comput.

// ==================
// Initializers
// ==================
ThirdOrderGauss::ThirdOrderGauss(int n, VectorXd initialMean, MatrixXd initialVariance,
                                 std::vector<double> param,
                                 VectorXd (*drift) (VectorXd, std::vector<double>&),
                                 MatrixXd (*diffusion) (VectorXd, std::vector<double>&, double, double),
                                 double dataVariance, double step, double sigmadiff, bool isStiff,
                                 double stiffIndex, double dampingInput)
{
    size = n;
    m = initialMean;
    oldM = m;
    P = initialVariance;
    f = drift;
    L = diffusion;
    xi = computeSigmaPoints();
    theta = param;
    varData = dataVariance * MatrixXd::Identity(size, size);
    h = step;
    sigma = sigmadiff;

    // Stabilized method
    stableMethod = isStiff;
    rho = stiffIndex;
    damping = dampingInput;
}

std::vector<VectorXd> ThirdOrderGauss::computeSigmaPoints(void)
{
    int n = size;
    std::vector<VectorXd> sigmaPoints(2 * n, VectorXd(n));
    double sqrtN = sqrt(n);

    for (int i = 0; i < n; i++) {
        sigmaPoints[i] = VectorXd::Zero(n);
        sigmaPoints[i + n] = VectorXd::Zero(n);
        sigmaPoints[i](i) = sqrtN;
        sigmaPoints[i + n](i) = -sqrtN;
    }

    return sigmaPoints;
}

// ==================
// Tools
// ==================
void ThirdOrderGauss::updateSqrtP(void)
{
    if (P == MatrixXd::Zero(size, size)) {
        sqrtP = MatrixXd::Zero(size, size);
    } else {
        sqrtP = P.sqrt();
    }
}

MatrixXd ThirdOrderGauss::matrixSqrt(MatrixXd& A)
{
    if (A == MatrixXd::Zero(size, size)) {
        return MatrixXd::Zero(size, size);
    } else {
        return A.sqrt();
    }
}

double ThirdOrderGauss::evaluateGaussian(VectorXd data)
{
    // SINCE WE ARE COMPUTING RATIOS, DO NOT INCLUDE THE MULTIPLICATIVE TERM
    double A = -0.5 * (data - m).transpose() * (varData + P).inverse() * (data - m);
    // std::cout << m.transpose() << std::endl << data.transpose() << std::endl << std::endl;
    return exp(A);
}

// ==================
// mean EFupdates
// ==================
VectorXd ThirdOrderGauss::meanUpdateFct(VectorXd mIntern, MatrixXd sqrtPIntern)
{
    VectorXd sum = VectorXd::Zero(size);

    double weight = 1.0 / (2 * size);

    for (int i = 0; i < 2 * size; i++) {
        sum +=  f(mIntern + sqrtPIntern * xi[i], theta);
    }
    sum *= weight;

    return sum;
}

// ==================
// variance EFupdates
// ==================
MatrixXd ThirdOrderGauss::varianceUpdateFct(VectorXd mIntern, MatrixXd sqrtPIntern)
{
    MatrixXd sumOne = MatrixXd::Zero(size, size);
    MatrixXd sumTwo = MatrixXd::Zero(size, size);
    MatrixXd sumThree = MatrixXd::Zero(size, size);

    double weight = 1.0 / (2 * size);

    VectorXd F(size);
    MatrixXd l(size, size);

    for (int i = 0; i < 2 * size; i++) {
        F = f(mIntern + sqrtPIntern * xi[i], theta);
        sumOne += F * xi[i].transpose() * sqrtPIntern.transpose();
        sumTwo += sqrtPIntern * xi[i] * F.transpose();
        l = L(mIntern + sqrtPIntern * xi[i], theta, sigma, h);
        sumThree += l * l.transpose();
    }

    return (sumOne + sumTwo + sumThree) * weight;
}

// ==================
// Un-stiff time integration
// Implemented with Euler Forward.
// ==================
void ThirdOrderGauss::EFupdates(int nSteps)
{
    updateSqrtP();
    for (int i = 0; i < nSteps; i++) {
        oldM = m;
        m += meanUpdateFct(m, sqrtP) * h;
        P += varianceUpdateFct(oldM, sqrtP) * h;
        updateSqrtP();
    }
}

// TODO: Debug this
void ThirdOrderGauss::stabUpdates(int nSteps)
{
    MatrixXd sqrtK(size, size);
    for (int j = 0; j < nSteps; j++) {
        // Hard-code the first two stages
        kMean[0] = m * stageCoeff[0][0];
        kVar[0] = P * stageCoeff[0][0];
        sqrtK = matrixSqrt(kVar[0]);
        kMean[1] = m * stageCoeff[1][0] + meanUpdateFct(kMean[0], sqrtK) * stageCoeff[1][1];
        kVar[1] = P * stageCoeff[1][0] + varianceUpdateFct(kMean[0], sqrtK) * stageCoeff[1][1];
        sqrtK = matrixSqrt(kVar[1]);
        for (int i = 2; i < nStages + 1; i++) {
            kMean[i] = meanUpdateFct(kMean[i - 1], sqrtK) * stageCoeff[i][0] +
                       kMean[i - 1] * stageCoeff[i][1] +
                       kMean[i - 2] * stageCoeff[i][2];
            kVar[i] = varianceUpdateFct(kMean[i - 1], sqrtK) * stageCoeff[i][0] +
                      kVar[i - 1] * stageCoeff[i][1] +
                      kVar[i - 2] * stageCoeff[i][2];
            sqrtK = matrixSqrt(kVar[i]);
        }
        m = kMean.back();
        P = kVar.back();
    }
}

// ==================
// Kalman updates
// ==================
void ThirdOrderGauss::KalmanUpdate(VectorXd data)
{
    MatrixXd S = P + varData;
    MatrixXd K = P * S.inverse();
    m = m + K * (data - m);
    P = P - K * S * K.transpose();
}

// ==================
// Final algorithm
// ==================
double ThirdOrderGauss::oneStep(VectorXd data, int nSteps)
{
    if (stableMethod) {
        nStages = static_cast<unsigned int> (std::max(std::ceil(sqrt(3.0 * rho / nSteps)), 2.0));
        std::cout << nStages << std::endl;
        omegaZero = 1.0 + damping / (nStages * nStages);
        chebCoeff.resize(nStages + 1);
        computeChebCoeff(omegaZero);
        stageCoeff.resize(nStages + 1);
        computeStageCoeff();
        kMean.resize(nStages + 1);
        kVar.resize(nStages + 1);
        for (int i = 0; i < nStages + 1; i++) {
            kMean[i].resize(size);
            kVar[i].resize(size, size);
        }
        stabUpdates(nSteps);
    } else {
        EFupdates(nSteps);
    }
    double result = evaluateGaussian(data);
    KalmanUpdate(data);
    return result;
}

// ================
// Stability update tools
// ================
void ThirdOrderGauss::computeChebCoeff(double x)
{
    chebCoeff[0] = 1.0;
    chebCoeff[1] = x;
    for (int i = 2; i < nStages + 1; i++) {
        chebCoeff[i] = 2 * x * chebCoeff[i - 1] - chebCoeff[i - 2];
    }

    // Compute omegaOne (derivative of Cheb polynomials)
    std::vector<double> chebDerivatives(nStages + 1);
    chebDerivatives[0] = 0.0;
    chebDerivatives[1] = 1.0;
    chebDerivatives[2] = 4 * x;
    double doubI;
    for (int i = 3; i < nStages + 1; i++) {
        doubI = static_cast<double>(i);
        chebDerivatives[i] = doubI * (2 * chebCoeff[i - 1] + chebDerivatives[i - 2] / (doubI - 2));
    }
    omegaOne = chebCoeff.back() / chebDerivatives.back();
}

void ThirdOrderGauss::computeStageCoeff(void)
{
    // Hard-code the first two stages
    stageCoeff[0].resize(1);
    stageCoeff[0][0] = 1.0;

    stageCoeff[1].resize(2);
    stageCoeff[1][0] = 1.0;
    stageCoeff[1][1] = h * omegaOne / omegaZero;

    for (int i = 2; i < nStages + 1; i++) {
        stageCoeff[i].resize(3);
        stageCoeff[i][0] = 2 * h * omegaOne * chebCoeff[i - 1] / chebCoeff[i];
        stageCoeff[i][1] = 2 * omegaZero * chebCoeff[i - 1] / chebCoeff[i];
        stageCoeff[i][2] = -1.0 * chebCoeff[i - 2] / chebCoeff[i];
    }
}

