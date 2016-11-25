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
    Fx = odeModel.odeJac;
    ODE = odeModel;
    L = diffusion;
    oldM = m;
    P = initialVariance;
    LLT<MatrixXd> chol(P);
    lChol = chol.matrixL();
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
        LLT<MatrixXd> chol(P);
        sqrtP = chol.matrixL();
    }
}

MatrixXd ThirdOrderGauss::matrixSqrt(MatrixXd& A)
{
    if (A == MatrixXd::Zero(size, size)) {
        return MatrixXd::Zero(size, size);
    } else {
        LLT<MatrixXd> chol(A);
        return chol.matrixL();
    }
}

double ThirdOrderGauss::evaluateGaussian(VectorXd data)
{
    /* std::cout << m.transpose() << std::endl << data.transpose() << std::endl
              << std::endl << P << std::endl << "================" << std::endl; */
    double A = -0.5 * (data - m).transpose() * (varData + P).inverse() * (data - m);
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

    double weight = 1.0 / (2 * size);

    VectorXd F(size);
    MatrixXd l(size, size);

    for (int i = 0; i < 2 * size; i++) {
        F = f(mIntern + sqrtPIntern * xi[i], theta);
        sumOne += F * xi[i].transpose() * sqrtPIntern.transpose();
        sumTwo += sqrtPIntern * xi[i] * F.transpose();
    }

    l = L(mIntern, theta, sigma, h);
    return (sumOne + sumTwo) * weight + l * l.transpose();
}

MatrixXd ThirdOrderGauss::TvarianceUpdateFct(VectorXd mIntern, MatrixXd PIntern)
{
    MatrixXd FxEval = Fx(mIntern, theta);
    MatrixXd l = L(mIntern, theta, sigma, h);
    return PIntern * FxEval.transpose() + FxEval * PIntern + l * l.transpose();
}

MatrixXd ThirdOrderGauss::TsqrtVarianceUpdateFct(VectorXd mIntern, MatrixXd W)
{
    MatrixXd FxEval = Fx(mIntern, theta);
    MatrixXd l = L(mIntern, theta, sigma, h);
    MatrixXd S = l * l.transpose();
    MatrixXd Sbar = S * 0.5;

    // Construction of an antisymmetric X as in Andrews (1968)
    /* MatrixXd LInv = triInv(LIntern, size);

    double sumOne, sumTwo, sumThree, sumFour;
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
    xAntiSym = xAntiSym - tmp; */

    // Square root update as in Tapley - Choe(1968)
    MatrixXd T = MatrixXd::Zero(size, size);
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < i; j++) {
            for (int k = j; k < size; k++) {
                T(i, j) += FxEval(i, k) * W(k, j);
            }
        }
    }

    MatrixXd M = MatrixXd::Zero(size, size);
    for (int i = 0; i < size; i++) {
        M(i, i) = Sbar(i, i);
        for (int k = 0; k < i - 1; k++) {
            M(i, i) -= M(i, k) * W(i, k);
        }
        M(i, i) = M(i, i) / W(i, i);
        for (int j = 0; j < i; j++) {
            for (int k = 0; k < j; k++) {
                M(i, j) -= M(j, k) * W(i, k);
            }
            for (int k = j; k < i; k++) {
                M(i, j) += T(j, k) * W(i, k);
            }
            for (int k = 0; k < j - 1; k++) {
                M(i, j) -= M(i, k) * W(j, k);
            }
            M(i, j) = M(i, j) / W(j, j);
        }
    }

    return MatrixXd((M + T).triangularView<Lower>());

    // Just choose C = 0
    /* MatrixXd triangularInverse = W.inverse();
    return FxEval * W + Sbar * triangularInverse.transpose(); */

    // Debug
    /* std::cout << LIntern
              << std::endl
              << "========"
              << std::endl
              << Sbar
              << std::endl
              << "========"
              << std::endl
              << M
              << std::endl
              << "========"
              << std::endl
              << T
              << std::endl
              << "======================"
              << std::endl; */
}

// ==================
// Un-stiff time integration
// Implemented with Euler Forward.
// ==================
void ThirdOrderGauss::EFupdates(int nSteps)
{
    for (int i = 0; i < nSteps; i++) {
        oldM = m;
        m += TmeanUpdateFct(m) * h;
        P += TvarianceUpdateFct(oldM, P) * h;
    }
    std::cout << "END"
              << std::endl
              << m.transpose()
              << std::endl
              << std::endl
              << P
              << std::endl
              << std::endl;
}

void ThirdOrderGauss::stabUpdates(int nSteps)
{
    if (stabParam.method == ROCK) {

        /* for (int j = 0; j < nSteps; j++) {
            std::cout << j * h << std::endl;
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
        } */

    } else if (stabParam.method == stdRKC) {
        for (int j = 0; j < nSteps; j++) {

            // Compute the stiffness index
            if (j % 1 == 0) {

                // For vdPol we know an analytic expression of the e'values...
                if (ODE.ode == VDPOL) {
                    double tmp = 1 - m(0) * m(0);
                    double delta = theta[0] * theta[0] * tmp * tmp -
                                   4 * theta[0] * (2 * m(0) * m(1) - 1);
                    double re = theta[0] * tmp;
                    if (delta > 0) {
                        double sqrtD = sqrt(delta);
                        double lambdaOne = (re + sqrtD) / 2;
                        double lambdaTwo = (re - sqrtD) / 2;
                        stabParam.stiffIndex = 2.0 * std::abs(std::min(lambdaOne, lambdaTwo));
                    } else {
                        stabParam.stiffIndex = 2.0 * std::abs(re / 2);
                    }
                } else {
                    // for the other equations, just use a power method
                    stabParam.stiffIndex = 4.0 * std::abs(powerMethod(Fx(m, theta), 1.0, 100));
                }
                stabParam.nStages = 2 + static_cast<int>(std::max(std::ceil(sqrt(0.5 * stabParam.stiffIndex * h)),
                                                                  2.0));

                /* std::cout << "time = "
                          << j * h
                          << " solution = "
                          << m.transpose()
                          << std::endl
                          << "stiffness = "
                          << stabParam.stiffIndex
                          << " nStagesRKC = "
                          << stabParam.nStages
                          << std::endl; */

                kMean.resize(stabParam.nStages + 1);
                kVar.resize(stabParam.nStages + 1);

                for (int i = 0; i < stabParam.nStages + 1; i++) {
                    kMean[i].resize(size);
                    kVar[i].resize(size, size);
                }
            }

            kMean[0] = m;
            kVar[0] = P;
            double coeff = h / static_cast<double>(stabParam.nStages * stabParam.nStages);
            kMean[1] = m + TmeanUpdateFct(m) * coeff;
            kVar[1] = P + TvarianceUpdateFct(m, P) * coeff;
            coeff = 2.0 * coeff;

            for (int i = 2; i < stabParam.nStages + 1; i++) {
                kMean[i] = TmeanUpdateFct(kMean[i - 1]) * coeff +
                           kMean[i - 1] * 2.0 - kMean[i - 2];
                kVar[i] = TvarianceUpdateFct(kMean[i - 1], kVar[i - 1]) * coeff +
                          kVar[i - 1] * 2.0 - kVar[i - 2];
            }

            m = kMean.back();
            P = kVar.back();


            /* kMean[0] = m;
            kVar[0] = lChol;
            double coeff = h / static_cast<double>(stabParam.nStages * stabParam.nStages);
            kMean[1] = m + TmeanUpdateFct(m) * coeff;
            kVar[1] = lChol + TsqrtVarianceUpdateFct(m, lChol) * coeff;
            coeff = 2.0 * coeff;

            for (int i = 2; i < stabParam.nStages + 1; i++) {
                kMean[i] = TmeanUpdateFct(kMean[i - 1]) * coeff +
                           kMean[i - 1] * 2.0 - kMean[i - 2];
                kVar[i] = TsqrtVarianceUpdateFct(kMean[i - 1], kVar[i - 1]) * coeff +
                          kVar[i - 1] * 2.0 - kVar[i - 2];
            }
            m = kMean.back();
            lChol = kVar.back(); */

            std::cout << m
                      << std::endl
                      << std::endl
                      << P
                      << std::endl
                      << "========="
                      << std::endl;

        }
        // P = lChol * lChol.transpose();

        std::cout << "END"
                  << std::endl
                  << m.transpose()
                  << std::endl
                  << std::endl
                  << P
                  << std::endl
                  << std::endl;
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
    LLT<MatrixXd> chol(P);
    lChol = chol.matrixL();
}

// ==================
// Final algorithm
// ==================
double ThirdOrderGauss::oneStep(VectorXd data, int nSteps)
{
    if (stableMethod) {
        if (stabParam.method == ROCK) {
            /* stabParam.nStages = static_cast<int> (std::max(std::ceil(sqrt(3.0 * stabParam.stiffIndex / nSteps)), 2.0));
            omegaZero = 1.0 + stabParam.damping / static_cast<double>(stabParam.nStages * stabParam.nStages);
            chebCoeff.resize(stabParam.nStages + 1);
            computeChebCoeff(omegaZero);
            stageCoeff.resize(stabParam.nStages + 1);
            computeStageCoeff();
            kMean.resize(stabParam.nStages + 1);
            kVar.resize(stabParam.nStages + 1);
            for (int i = 0; i < stabParam.nStages + 1; i++) {
                kMean[i].resize(size);
                kVar[i].resize(size, size); */
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

// Get functions
VectorXd ThirdOrderGauss::getMean(void)
{
    return m;
}

MatrixXd ThirdOrderGauss::getVariance(void)
{
    return P;
}