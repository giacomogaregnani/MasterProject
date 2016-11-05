#include <cmath>
#include <random>
#include <iostream>
#include "mcmcTools.hpp"
#include <Solver.hpp>

#define PI 3.1415926535897	
#define M_EPS 0.0000001

// MLCMWMCMC
std::vector<VectorXd> sMLmetropolisHastings(VectorXd initialCond, std::vector<double>& param, double sigma, 
					    int size, double h, double finalTime, 
					    std::vector<VectorXd>& data, 
					    std::vector<double>& dataTimes, 
					    VectorXd (*odeFunc) (VectorXd, std::vector<double>&), 
					    VectorXd priorMean, VectorXd priorVariance, int internalMC, int nStepsMC, double varData,
				            double rho, double damping)
{
	double gamma = 0.01;
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
	double l;
	double prior;
	double oldParamL;
	double oldLike;
	double oldPrior = evalLogPrior(paramVec, priorMean, priorVariance, nParam); 
	
	oldLike = 0.0;
	for (size_t i = 0; i < nData; i++) {
		sMLMC<detSROCK> multilevel(size, h, initialCond, param, odeFunc, sigma, 1, dataTimes[i], 
			                   damping, rho, true, h, &evalSingleLikelihood, data[i], varData);
		oldLike += multilevel.compute();	
	}
	
	VectorXd paramVecN(nParam); 	
	std::vector<VectorXd> mcmcPath(nStepsMC, VectorXd(nParam));
	VectorXd w(nParam);
	std::vector<double> paramStd(nParam);
	
	for (int i = 0; i < nStepsMC; i++) {
		if (i % 50 == 0) {
			std::cout << "iteration: " << i << std::endl;
			std::cout << "likelihood: " << oldLike << std::endl;
			std::cout << "step: " << gamma << std::endl;
			std::cout << "prior: " << oldPrior << std::endl;
			std::cout << "param: " << paramVec.transpose() << std::endl;
		}

		// Generate new guess for parameter
		for (int j = 0; j < nParam; j++) {
			w(j) = gamma * normal(wGenerator);	
		}
		paramVecN = paramVec + w;
		for (int j = 0; j < nParam; j++) {
			paramStd[j] = paramVecN(j);
		}

		l = 0.0;	
		// MLMC for the new value of the parameter
		for (size_t j = 0; j < nData; j++) {
			sMLMC<detSROCK> multilevel(size, h, initialCond, paramStd, odeFunc, sigma, 1, dataTimes[j], 
				                   damping, rho, true, h, &evalSingleLikelihood, data[j], varData);
			l += multilevel.compute();
		}
		prior = evalLogPrior(paramVecN, priorMean, priorVariance, nParam);
		
		// MLMC for the old value of the parameter
		for (int j = 0; j < nParam; j++) {
			paramStd[j] = paramVec(j);
		}
	
		oldParamL = 0.0;	
		for (size_t j = 0; j < nData; j++) {
			sMLMC<detSROCK> multilevel(size, h, initialCond, paramStd, odeFunc, sigma, 1, dataTimes[j], 
				                   damping, rho, true, h, &evalSingleLikelihood, data[j], varData);
			oldParamL += multilevel.compute();	
		}

		// Generate probability and update 
		double u = log(unif(wGenerator));
		double alpha = std::min<double>(0.0, (l + prior) - (oldParamL + oldPrior));	
		if (u < alpha) {
			paramVec = paramVecN;
			oldLike = l;
			oldPrior = prior;
			gamma = std::min<double>(gamma * 2.0, 1.0);
		} else {
			gamma = std::max<double>(gamma / 2.0, 0.01);
		}	
		mcmcPath[i] = paramVec;	
	}
	return mcmcPath;
}
