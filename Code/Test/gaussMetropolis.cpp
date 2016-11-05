#include <cmath>
#include <random>
#include <iostream>
#include "mcmcTools.hpp"
#include <Solver.hpp>
#include <Eigen/Cholesky>

#define PI 3.1415926535897	
#define M_EPS 0.0000001

MatrixXd diffusion(VectorXd x, std::vector<double>& param, double sigma, double h)
{
	int size = x.size();
	MatrixXd M = sigma * h * MatrixXd::Identity(size, size); 
	return M;
}

std::vector<VectorXd> gaussMetropolisHastings(VectorXd initialCond, std::vector<double>& param, double sigma, 
			    int size, double h, double finalTime, std::vector<VectorXd>& data,
			    std::vector<double>& dataTimes, 
			    VectorXd (*odeFunc) (VectorXd, std::vector<double>&), 
			    VectorXd priorMean, VectorXd priorVariance, int nStepsMC, double varData, double* accRatio)
{
	*accRatio = 0.0;
	std::default_random_engine wGenerator{(unsigned int) time(NULL)};
	std::default_random_engine solGenerator{(unsigned int) time(NULL)};
	std::normal_distribution<double> normal(0.0, 1.0);
	std::uniform_real_distribution<double> unif;
	size_t nData = data.size();
	
	int nParam = (int) param.size();
	VectorXd paramVec(nParam);
	for (int i = 0; i < nParam; i++) {
		paramVec(i) = param[i];
	} 
	double l, prior, oldParamL, oldLike;

	// compute number of step for going from data to data with time step h
	std::vector<int> nSteps(nData);
	nSteps[0] = static_cast<int>(round(dataTimes[0] / h));
	for (size_t i = 1; i < nData; i++) {
		nSteps[i] = static_cast<int>(round((dataTimes[i] - dataTimes[i - 1]) / h));
	}

	oldLike = 1.0;
	MatrixXd initialVar = MatrixXd::Zero(size, size);
	ThirdOrderGauss solver(size, initialCond, initialVar, param, odeFunc, &diffusion, varData, h, sigma, true, 30.0, 0.5);
	for (size_t i = 0; i < nData; i++) {
		oldLike *= solver.oneStep(data[i], nSteps[i]);	
	}
	double oldPrior = evalPrior(paramVec, priorMean, priorVariance, nParam);


	VectorXd paramVecN(nParam); 	
	std::vector<VectorXd> mcmcPath(nStepsMC, VectorXd(nParam));
	VectorXd w(nParam);
	std::vector<double> paramStd(nParam);

	// =========== RAM UPDATE ===========
	double gamma = 0.01;
	MatrixXd I = MatrixXd::Identity(nParam, nParam);
	MatrixXd initialCov = gamma * I;
	LLT<MatrixXd> chol(initialCov);
	MatrixXd S = chol.matrixL();
	double desiredAlpha = 0.25;
	MatrixXd tmp(nParam, nParam);
	// ==================================

	for (int i = 0; i < nStepsMC; i++) {
		if (i % 50 == 0) {
			std::cout << "iteration: " << i << std::endl;
			std::cout << "likelihood: " << oldLike << std::endl;
			std::cout << "S matrix: " << std::endl << S << std::endl;
			std::cout << "prior: " << oldPrior << std::endl;
			std::cout << "param: " << paramVec.transpose() << std::endl;
		}

        // std::cout << "iteration " << i << " parameter " << paramVec(0) << std::endl;

		// Generate new guess for parameter
		for (int j = 0; j < nParam; j++) {
			w(j) = normal(wGenerator);
		}
		paramVecN = paramVec + S * w;

		// Gauss for the new value of the parameter
		for (int j = 0; j < nParam; j++) {
			paramStd[j] = paramVecN(j);
		}

		l = 1.0;
		ThirdOrderGauss newSolver(size, initialCond, initialVar, paramStd, odeFunc, &diffusion, varData, h, sigma, true, 30.0, 0.5);
		for (size_t j = 0; j < nData; j++) {
			l *= newSolver.oneStep(data[j], nSteps[j]);
		}
		prior = evalPrior(paramVecN, priorMean, priorVariance, nParam);
		
		// Gauss for the old value of the parameter
		for (int j = 0; j < nParam; j++) {
			paramStd[j] = paramVec(j);
		}
	
		oldParamL = 1.0;	
		ThirdOrderGauss oldSolver(size, initialCond, initialVar, paramStd, odeFunc, &diffusion, varData, h, sigma, true, 30.0, 0.5);
		for (size_t j = 0; j < nData; j++) {
			oldParamL *= oldSolver.oneStep(data[j], nSteps[j]);	
		}

		// Generate probability and update 
		double u = unif(wGenerator);
		double alpha = std::max<double>(0.0, std::min<double>(1.0, (l * prior) / (oldParamL * oldPrior)));
		if (u < alpha) {
			*accRatio += 1.0 / nStepsMC;
			paramVec = paramVecN;
			oldLike = l;
			oldPrior = prior;
			gamma = std::min<double>(gamma * 2.0, 1.0);
		} else {
			gamma = std::max<double>(gamma / 2.0, 0.01);
		}	
		mcmcPath[i] = paramVec;

		// RAM UPDATE
		MatrixXd WWT = w * w.transpose();
		double WTW = w.dot(w);
		double gammaI = std::min(1.0, 2.0 * pow(static_cast<double>(i + 1), -0.75));
		double diffAlpha = alpha - desiredAlpha;
		double coeff = gammaI * diffAlpha / WTW;
		MatrixXd C = I + coeff * WWT;
		tmp = S * C * S.transpose();
		LLT<MatrixXd> cholIt(tmp);
		S = cholIt.matrixL();
	}
	return mcmcPath;
}
