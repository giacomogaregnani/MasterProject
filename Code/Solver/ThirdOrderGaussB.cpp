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
                                 double dataVariance, double step, double sigmadiff)
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

double ThirdOrderGauss::evaluateGaussian(VectorXd data)
{
	// SINCE WE ARE COMPUTING RATIOS, DO NOT INCLUDE THE MULTIPLICATIVE TERM
	double A = -0.5 * (data - m).transpose() * (varData + P).inverse() * (data - m);
    // std::cout << m.transpose() << std::endl << data.transpose() << std::endl << std::endl;
	return exp(A);
}

// ==================
// mean updates 
// ==================
VectorXd ThirdOrderGauss::meanUpdateFct(void)
{
	VectorXd sum = VectorXd::Zero(size);
	
	double weight = 1.0 / (2 * size);

	for (int i = 0; i < 2 * size; i++) {
		sum +=  f(m + sqrtP * xi[i], theta); 
	}
	sum *= weight;

	return sum;		
}

// ==================
// variance updates 
// ==================
MatrixXd ThirdOrderGauss::varianceUpdateFct(void)
{
	MatrixXd sumOne = MatrixXd::Zero(size, size);
	MatrixXd sumTwo = MatrixXd::Zero(size, size);
	MatrixXd sumThree = MatrixXd::Zero(size, size);

	double weight = 1.0 / (2 * size);

	VectorXd F(size);
	MatrixXd l(size, size);

	for (int i = 0; i < 2 * size; i++) {
		F = f(oldM + sqrtP * xi[i], theta);	
		sumOne += F * xi[i].transpose() * sqrtP.transpose(); 
		sumTwo += sqrtP * xi[i] * F.transpose();
		l = L(oldM + sqrtP * xi[i], theta, sigma, h);		
		sumThree += l * l.transpose();
	}

	return weight * (sumOne + sumTwo + sumThree);
}

// ==================
// time integration
// Implemented with Euler Forward.
// TODO: Explore other methods
// ==================
void ThirdOrderGauss::updates(int nSteps)
{
	updateSqrtP();
	for (int i = 0; i < nSteps; i++) {
		oldM = m;
        m += h * meanUpdateFct();
		P += h * varianceUpdateFct();
		updateSqrtP();
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
	updates(nSteps);	
	double result = evaluateGaussian(data);
	KalmanUpdate(data);
	return result;
}
