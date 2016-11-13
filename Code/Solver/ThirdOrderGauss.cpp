#include "Solver.hpp"
#include <unsupported/Eigen/MatrixFunctions>
#include <Eigen/Cholesky>
#include <Eigen/LU>
#include <cmath>
#include <fstream>
#include <string>

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
ThirdOrderGauss::ThirdOrderGauss(odeDef odeModel, MatrixXd initialVariance,
                                 std::vector<double> param,
                                 MatrixXd (*diffusion) (VectorXd, std::vector<double>&, double, double),
                                 double dataVariance, double step, double sigmadiff, bool isStiff,
                                 StabValues* stability)
{
    size = odeModel.size;
    f = odeModel.odeFunc;
    m = odeModel.initialCond;
    oldM = m;
    P = initialVariance;
    LLT<MatrixXd> chol(P);
    lChol = chol.matrixL();
    L = diffusion;
    xi = computeSigmaPoints();
    theta = param;
    varData = MatrixXd::Identity(size, size) * dataVariance;
    h = step;
    sigma = sigmadiff;
    xAntiSym = MatrixXd::Zero(size, size);

    // Stabilized method
    stableMethod = isStiff;
    if (stableMethod) {
        stabParam = *stability;
    }
    Fx = odeModel.odeJac;
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
    /* std::cout << data.transpose()
              << std::endl
              << m.transpose()
              << std::endl
              << P
              << std::endl
              << theta[0]
              << std::endl
              << "============"
              << std::endl; */
    if (A > 0) {
        throw 1;
    }
    return A;
}

// ==================
// mean updates
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

VectorXd ThirdOrderGauss::TmeanUpdateFct(VectorXd mIntern)
{
    return f(mIntern, theta);
}

// ==================
// variance updates
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

MatrixXd ThirdOrderGauss::TvarianceUpdateFct(VectorXd mIntern, MatrixXd PIntern)
{
    MatrixXd FxEval = Fx(mIntern, theta);
    MatrixXd l = L(mIntern, theta, sigma, h);
    return PIntern * FxEval.transpose() + FxEval * PIntern + l * l.transpose();
}

MatrixXd ThirdOrderGauss::TsqrtVarianceUpdateFct(VectorXd mIntern, MatrixXd LIntern)
{
    MatrixXd FxEval = Fx(mIntern, theta);
    MatrixXd l = L(mIntern, theta, sigma, h);
    MatrixXd S = l * l.transpose();
    MatrixXd LInv = triInv(LIntern, size);

    // Construction of an antisymmetric X as in Andrews (1968)
    /* double sumOne, sumTwo, sumThree, sumFour;
    for (int jX = 0; jX < size; jX++) {
        for (int iX = 0; iX < jX; iX++) {
            sumOne = 0.0;
            for (int kX = jX; kX < size; kX++) {
                sumOne += FxEval(iX, kX) * LIntern(kX, jX);
            }
            sumTwo = 0.0;
            for (int mX = 0; mX < jX + 1; mX++) {
                sumTwo += S(iX, mX) * LInv(jX, mX);
            }
            sumTwo *= 0.5;
            sumThree = 0.0;
            for (int mX = 0; mX < iX; mX++) {
                sumThree += xAntiSym(iX, mX) * LInv(jX, mX);
            }
            sumFour = 0.0;
            for (int mX = iX + 1; mX < jX; mX++) {
                sumFour += xAntiSym(mX, iX) * LInv(jX, mX);
            }
            xAntiSym(jX, iX) = LIntern(jX, jX) * (sumOne + sumTwo + sumThree - sumFour);
        }
    }
    MatrixXd tmp = xAntiSym.transpose();
    xAntiSym = xAntiSym - tmp;

    std::cout << LIntern
              << std::endl
              << "======="
              << std::endl
              << LInv
              << std::endl
              << "======="
              << std::endl
              << xAntiSym
              << std::endl
              << "======="
              << std::endl
              << FxEval * LIntern + (xAntiSym + S * 0.5) * LInv.transpose()
              << std::endl
              << "======================"
              << std::endl;
    */

    return FxEval * LIntern + (xAntiSym + S * 0.5) * LInv.transpose();
    // return MatrixXd((FxEval * LIntern + (xAntiSym + S * 0.5) * LInv.transpose()).triangularView<Lower>());
}

// ==================
// Un-stiff time integration
// Implemented with Euler Forward.
// ==================
void ThirdOrderGauss::EFupdates(int nSteps)
{
    //updateSqrtP();
    /*for (int i = 0; i < nSteps; i++) {
        oldM = m;
        m += TmeanUpdateFct(m) * h;
        P += TvarianceUpdateFct(oldM, P) * h;
        //updateSqrtP();
    } */

    for (int i = 0; i < nSteps; i++) {
        oldM = m;
        m += TmeanUpdateFct(m) * h;
        lChol += TsqrtVarianceUpdateFct(oldM, lChol) * h;
    }
    P = lChol * lChol.transpose();
}

void ThirdOrderGauss::stabUpdates(int nSteps)
{
    // MatrixXd sqrtK(size, size);
    if (stabParam.method == ROCK) {
        for (int j = 0; j < nSteps; j++) {

            // SQRT ROCK TAYLOR UPDATE
            kMean[0] = m * stageCoeff[0][0];
            kVar[0] = lChol * stageCoeff[0][0];

            kMean[1] = m * stageCoeff[1][0] + TmeanUpdateFct(kMean[0]) * stageCoeff[1][1];
            kVar[1] = lChol * stageCoeff[1][0] + TsqrtVarianceUpdateFct(kMean[0], lChol) * stageCoeff[1][1];

            for (int i = 2; i < stabParam.nStages + 1; i++) {
                kMean[i] = TmeanUpdateFct(kMean[i - 1]) * stageCoeff[i][0] +
                           kMean[i - 1] * stageCoeff[i][1] +
                           kMean[i - 2] * stageCoeff[i][2];
                kVar[i] = TsqrtVarianceUpdateFct(kMean[i - 1], kVar[i - 1]) * stageCoeff[i][0] +
                          kVar[i - 1] * stageCoeff[i][1] +
                          kVar[i - 2] * stageCoeff[i][2];
            }
            m = kMean.back();
            lChol = kVar.back();
        }
    } else if (stabParam.method == stdRKC) {
        for (int j = 0; j < nSteps; j++) {

            // RKC SQRT TAYLOR UPDATE
            kMean[0] = m;
            kVar[0] = lChol;
            double coeff = h / static_cast<double>(stabParam.nStages * stabParam.nStages);
            kMean[1] = m + TmeanUpdateFct(m) * coeff;
            kVar[1] = lChol + TsqrtVarianceUpdateFct(m, lChol) * coeff;
            P = kVar[1] * kVar[1].transpose();
            LLT<MatrixXd> chol(P);
            kVar[1] = chol.matrixL();
            coeff = 2.0 * coeff;

            for (int i = 2; i < stabParam.nStages + 1; i++) {
                kMean[i] = TmeanUpdateFct(kMean[i - 1]) * coeff +
                           kMean[i - 1] * 2.0 - kMean[i - 2];
                kVar[i] = TsqrtVarianceUpdateFct(kMean[i - 1], kVar[i - 1]) * coeff +
                          kVar[i - 1] * 2.0 - kVar[i - 2];
                P = kVar[i] * kVar[i].transpose();
                LLT<MatrixXd> chol(P);
                kVar[i] = chol.matrixL();
            }

            m = kMean.back();
            lChol = kVar.back();

        }
    }
    P = lChol * lChol.transpose();
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
        if (stabParam.method == stdRKC) {
            stabParam.nStages = static_cast<int>(std::max(std::ceil(sqrt(0.5 * stabParam.stiffIndex / nSteps)), 2.0));
        } else if (stabParam.method == ROCK) {
            stabParam.nStages = static_cast<int> (std::max(std::ceil(sqrt(3.0 * stabParam.stiffIndex / nSteps)), 2.0));
            omegaZero = 1.0 + stabParam.damping / static_cast<double>(stabParam.nStages * stabParam.nStages);
            chebCoeff.resize(stabParam.nStages + 1);
            computeChebCoeff(omegaZero);
            stageCoeff.resize(stabParam.nStages + 1);
            computeStageCoeff();
        }
        kMean.resize(stabParam.nStages + 1);
        kVar.resize(stabParam.nStages + 1);
        for (int i = 0; i < stabParam.nStages + 1; i++) {
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
// Stability update tools (ROCK)
// ================
void ThirdOrderGauss::computeChebCoeff(double x)
{
    chebCoeff[0] = 1.0;
    chebCoeff[1] = x;
    for (int i = 2; i < stabParam.nStages + 1; i++) {
        chebCoeff[i] = 2.0 * x * chebCoeff[i - 1] - chebCoeff[i - 2];
    }

    // Compute omegaOne (derivative of Cheb polynomials)
    std::vector<double> chebDerivatives(stabParam.nStages + 1);
    chebDerivatives[0] = 0.0;
    chebDerivatives[1] = 1.0;
    chebDerivatives[2] = 4.0 * x;
    double doubI;
    for (int i = 3; i < stabParam.nStages + 1; i++) {
        doubI = static_cast<double>(i);
        chebDerivatives[i] = doubI * (2 * chebCoeff[i - 1] + chebDerivatives[i - 2] / (doubI - 2.0));
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

    for (int i = 2; i < stabParam.nStages + 1; i++) {
        stageCoeff[i].resize(3);
        stageCoeff[i][0] = 2 * h * omegaOne * chebCoeff[i - 1] / chebCoeff[i];
        stageCoeff[i][1] = 2 * omegaZero * chebCoeff[i - 1] / chebCoeff[i];
        stageCoeff[i][2] = -1.0 * chebCoeff[i - 2] / chebCoeff[i];
    }
}