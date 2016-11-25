#include <cmath>
#include <random>
#include <iostream>
#include "mcmcTools.hpp"
#include <Solver.hpp>

#define PI 3.1415926535897	

// MLCMWMCMC
std::vector<VectorXd> sMLmetropolisHastings(odeDef odeModel, std::vector<double>& param, double sigma,
											double h, double finalTime,
											std::vector<VectorXd>& data,
											std::vector<double>& dataTimes,
											VectorXd priorMean, VectorXd priorVariance,
											int nStepsMC, double varData, double damping, long int* cost)
{
	// Initialize cost to zero
	*cost = 0;

	// Create generators
	std::default_random_engine wGenerator{(unsigned int) time(NULL)};
	std::default_random_engine solGenerator{(unsigned int) time(NULL)};
	std::normal_distribution<double> normal(0.0, 1.0);
	std::uniform_real_distribution<double> unif;
	size_t nData = data.size();

	// Get data from odeModel
	int size = odeModel.size;
	VectorXd (*odeFunc) (VectorXd, std::vector<double>&) = odeModel.odeFunc;
	VectorXd initialCond = odeModel.initialCond;

	// Initialize parameter vector
	int nParam = (int) param.size();
	VectorXd paramVec(nParam);
	for (int i = 0; i < nParam; i++) {
		paramVec(i) = param[i];
	}

	// Initialize support doubles
	double l, prior, oldParamL, oldLike, oldGauss, newGauss;
	double oldPrior = evalLogPrior(paramVec, priorMean, priorVariance, nParam);

	// Compute stiffness index for initial guess of theta
	double rho;
	rho = 2.0 * std::abs(powerMethod(odeModel.odeJac(odeModel.initialCond, param), 1.0, 100));

	// Compute likelihood of inital guess
	oldLike = 0.0;
	for (size_t i = 0; i < nData; i++) {
		sMLMC<detSROCK> multilevel(size, h, initialCond, param, odeFunc, sigma, 1, dataTimes[i],
							  damping, rho, true, h, &evalSingleLikelihood, data[i], varData);
		oldLike += multilevel.compute();
		*cost += multilevel.cost();
	}

	// Initialize support structures
	VectorXd paramVecN(nParam); 	
	std::vector<VectorXd> mcmcPath(nStepsMC, VectorXd(nParam));
	VectorXd w(nParam);
	std::vector<double> paramStd(nParam);

	// Initialize RAM
	double gamma = 0.01;
	double desiredAlpha = 0.25;
	MatrixXd S = RAMinit(gamma, desiredAlpha, nParam);

	for (int i = 0; i < nStepsMC; i++) {
		if (i % 50 == 0) {
			std::cout << "iteration: " << i << std::endl;
			std::cout << "likelihood: " << oldLike << std::endl;
			std::cout << "step: " << S << std::endl;
			std::cout << "prior: " << oldPrior << std::endl;
			std::cout << "param: " << paramVec.transpose() << std::endl;
		}

		// Generate new guess for parameter
		bool isPositive = false;
		while (!isPositive) {
			for (int j = 0; j < nParam; j++) {
				w(j) = normal(wGenerator);
			}
			paramVecN = paramVec + S * w;
			if (paramVecN(0) > 0) {
				isPositive = true;
				oldGauss = phi(paramVec(0) / (S(0, 0) * S(0, 0)));
				newGauss = phi(paramVecN(0) / (S(0, 0) * S(0, 0)));
			}
		}
		for (int j = 0; j < nParam; j++) {
			paramStd[j] = paramVecN(j);
		}

		// MLMC for the new value of the parameter
		rho = 2.0 * std::abs(powerMethod(odeModel.odeJac(odeModel.initialCond, paramStd), 1.0, 100));
		l = 0.0;
		for (size_t j = 0; j < nData; j++) {
			sMLMC<detSROCK> multilevel(size, h, initialCond, paramStd, odeFunc, sigma, 1, dataTimes[j],
								  damping, rho, true, h, &evalSingleLikelihood, data[j], varData);
			l += multilevel.compute();
			*cost += multilevel.cost();
		}
		prior = evalLogPrior(paramVecN, priorMean, priorVariance, nParam);

		// MLMC for the old value of the parameter
		for (int j = 0; j < nParam; j++) {
			paramStd[j] = paramVec(j);
		}
		rho = 2.0 * std::abs(powerMethod(odeModel.odeJac(odeModel.initialCond, paramStd), 1.0, 100));
		oldParamL = 0.0;	
		for (size_t j = 0; j < nData; j++) {
			sMLMC<detSROCK> multilevel(size, h, initialCond, paramStd, odeFunc, sigma, 1, dataTimes[j],
				                   damping, rho, true, h, &evalSingleLikelihood, data[j], varData);
			oldParamL += multilevel.compute();
			*cost += multilevel.cost();
		}

		// Generate probability and update 
		double u = log(unif(wGenerator));
		double alpha = std::min<double>(0.0, (l + prior + log(oldGauss)) - (oldParamL + oldPrior + log(newGauss)));
		if (u < alpha) {
			paramVec = paramVecN;
			oldLike = l;
			oldPrior = prior;
		}
		mcmcPath[i] = paramVec;

		S = RAMupdate(S, w, exp(alpha), desiredAlpha, nParam, i + 1);
	}

	return mcmcPath;
}