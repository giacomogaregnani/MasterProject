#include <cmath>
#include <random>
#include <iostream>
#include "mcmcTools.hpp"
#include <Solver.hpp>
#include<Eigen/Cholesky>

#define PI 3.1415926535897
#define M_EPS 0.0000001

// MLMCWM
std::vector<VectorXd> MLmetropolisHastings(VectorXd initialCond, std::vector<double>& param, double sigma,
										   int size, double h, double finalTime, std::vector<VectorXd>& data,
										   std::vector<double>& dataTimes,
										   VectorXd (*odeFunc) (VectorXd, std::vector<double>&),
										   VectorXd priorMean, VectorXd priorVariance, int internalMC,
										   int nStepsMC, double varData,  long int* costPerIteration,
                                           double* acceptRatio)
{
    *acceptRatio = 0.0;
	std::default_random_engine wGenerator{(unsigned int) time(NULL)};
	std::default_random_engine solGenerator{(unsigned int) time(NULL)};
	std::normal_distribution<double> normal(0.0, 1.0);
	std::uniform_real_distribution<double> unif;
	size_t nData = data.size();

    int nParam = static_cast<int>(param.size());
	VectorXd paramVec(nParam);
	for (int i = 0; i < nParam; i++) {
		paramVec(i) = param[i];
	}
	double l, prior, oldParamL, oldLike;

	oldLike = 0.0;
	for (size_t i = 0; i < nData; i++) {
		MLMC<EulerForward> multilevel(size, h, initialCond, param, odeFunc, sigma, 1,
									  dataTimes[i], true, h, &evalSingleLikelihood, data[i], varData);
        oldLike += multilevel.compute();
        *(costPerIteration) = multilevel.cost();
	}
    double oldPrior = evalPrior(paramVec, priorMean, priorVariance, nParam);

    VectorXd paramVecN(nParam);
    VectorXd w(nParam);
    std::vector<VectorXd> mcmcPath(nStepsMC, VectorXd(nParam));
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

		// Generate new guess for parameter
		for (int j = 0; j < nParam; j++) {
			w(j) = normal(wGenerator);
		}
		paramVecN = paramVec + S * w;

		// MLMC - NEW PARAMETER GUESS
		for (int j = 0; j < nParam; j++) {
			paramStd[j] = paramVecN(j);
		}

		l = 0.0;
		for (size_t j = 0; j < nData; j++) {
			MLMC<EulerForward> multilevel(size, h, initialCond, paramStd, odeFunc, sigma, 1,
										  dataTimes[j], true, h, &evalSingleLikelihood, data[j], varData);
			l += multilevel.compute();
		}
		prior = evalPrior(paramVecN, priorMean, priorVariance, nParam);

        // MLMC - OLD PARAMETER GUESS (PRIOR ALREADY COMPUTED)
		for (int j = 0; j < nParam; j++) {
			paramStd[j] = paramVec(j);
		}

		oldParamL = 0.0;
		for (size_t j = 0; j < nData; j++) {
			MLMC<EulerForward> multilevel(size, h, initialCond, paramStd, odeFunc, sigma, 1,
										  dataTimes[j], true, h, &evalSingleLikelihood, data[j], varData);
			oldParamL += multilevel.compute();
		}

		// Generate probability and update
		double u = unif(wGenerator);
		double alpha = std::max<double>(std::min<double>(1.0, (l * prior) / (oldParamL * oldPrior)), 0.0);
        if (u < alpha) {
            *(acceptRatio) += 1.0 / nStepsMC;
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
