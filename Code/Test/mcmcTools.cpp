#include <cmath>
#include <random>
#include <iostream>
#include "mcmcTools.hpp"
#include <Solver.hpp>

#define PI 3.1415926535897	

double evalLikelihood(std::vector<VectorXd> data,
					  std::vector<VectorXd> values,
					  int size, double variance)
{
	double likelihood = 0.0,  temp; 
	size_t n = data.size();
	
	for (size_t i = 0; i < n; i++) {
		for (int j = 0; j < size; j++) {
			temp = data[i](j) - values[i](j);
			likelihood += temp * temp;
		}			
	}
	return exp(-1.0 / (2.0 * variance) * likelihood);
}

double evalLogLikelihood(std::vector<VectorXd>& data,
						 std::vector<VectorXd>& values,
						 int size, double variance)
{
	double likelihood = 0.0; 
	size_t n = data.size();
	
	for (size_t i = 0; i < n; i++) {
		likelihood += (data[i] - values[i]).dot(data[i] - values[i]);	
	}	

	return -likelihood / (2.0 * variance);
}

double evalSingleLikelihood(VectorXd data, VectorXd value, double variance)
{
    double lik = (data - value).dot(data - value);
	return -lik / (2.0 * variance);
}


double evalPrior(VectorXd x, VectorXd mean, VectorXd variance, int size)
{
	// normal
	double prior = 1.0;
	for (int i = 0; i < size; i++) {
		prior *= exp(-(x(i) - mean(i)) * (x(i) - mean(i)) / (2 * variance(i)));
	}
	// lognormal
	/*double prior = 1.0;
	double logar;
	for (int i = 0; i < size; i++) {
		if (x(i) < 0) return 0.0;
		else logar = log(x(i));
		prior *= exp(-(logar - mean(i)) * (logar - mean(i)) / (2.0 * variance(i)));
	} */
	return prior;
}

double evalLogPrior(VectorXd x, VectorXd mean, VectorXd variance, int size)
{
	// normal 
	return -(x - mean).dot(x - mean) / (2 * variance(0));

	// log-normal distribution
	/* VectorXd logar(size);
	for (int i = 0; i < size; i++) {
		if (x(i) < 0) {
			return 0.0;
		} else {
			logar(i) = log(x(i));
		}
	}
	return -(logar - mean).dot(logar - mean) / (2 * variance(0));*/
}

MatrixXd RAMinit(double gamma, double desiredAlpha, int nParam)
{
    MatrixXd I = MatrixXd::Identity(nParam, nParam);
    MatrixXd initialCov = gamma * I;
    LLT<MatrixXd> chol(initialCov);
    return chol.matrixL();
}

MatrixXd RAMupdate(MatrixXd S, VectorXd w, double alpha, double desiredAlpha, int nParam, int nIt)
{
    MatrixXd WWT = w * w.transpose();
    double WTW = w.dot(w);
    double gammaI = std::min(1.0, 2.0 * pow(static_cast<double>(nIt), -0.75));
    double diffAlpha = alpha - desiredAlpha;
    double coeff = gammaI * diffAlpha / WTW;
    MatrixXd C = MatrixXd::Identity(nParam, nParam) + WWT * coeff;
    LLT<MatrixXd> cholIt(S * C * S.transpose());
    return cholIt.matrixL();
}

double phi(double x)
{
	// constants
	double a1 =  0.254829592;
	double a2 = -0.284496736;
	double a3 =  1.421413741;
	double a4 = -1.453152027;
	double a5 =  1.061405429;
	double p  =  0.3275911;

	// Save the sign of x
	int sign = 1;
	if (x < 0)
		sign = -1;
	x = fabs(x)/sqrt(2.0);

	// A&S formula 7.1.26
	double t = 1.0/(1.0 + p*x);
	double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);

	return 0.5*(1.0 + sign*y);
}