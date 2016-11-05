#include <cmath>
#include <random>
#include <iostream>
#include "mcmcTools.hpp"
#include <Solver.hpp>

#define PI 3.1415926535897	
#define M_EPS 0.0000001

std::vector<VectorXd> sMetropolisHastings(VectorXd initialCond, std::vector<double>& param, double sigma, 
			    int size, double h, double finalTime, 
                            std::vector<VectorXd>& data, 
			    std::vector<double>& dataTimes, 
			    VectorXd (*odeFunc) (VectorXd, std::vector<double>&), 
			    VectorXd priorMean, VectorXd priorVariance, int internalMC, 
			    int nStepsMC, int nLevels, double damping, double varData)
{
	double gamma = 0.01;
	std::default_random_engine wGenerator{(unsigned int) time(NULL)};
	std::default_random_engine solGenerator{(unsigned int) time(NULL)};
	std::normal_distribution<double> normal(0.0, 1.0);
	std::uniform_real_distribution<double> unif;
	size_t nData = data.size();
	
	double oldLike = 0.0;
	int nParam = (int) param.size();
	VectorXd paramVec(nParam);
	for (int i = 0; i < nParam; i++) {
		paramVec(i) = param[i];
	} 
	double l;
	double prior;

	// Evaluate posterior for initial guess
	sProbMethod<detSROCK> firstSolver(size, h, initialCond, param, odeFunc, sigma, nLevels, damping); 		
	for (int j = 0; j < internalMC; j++) {	
		double t = 0; 
		int count = 0;
		std::vector<VectorXd> solutionAtDataTimes(nData, VectorXd(size));
		while (t < finalTime) {
			firstSolver.oneStep(solGenerator, h);
			t = t + h;
			if (std::abs(t - dataTimes[count]) < h / 10.0) {
				std::cout << t << std::endl;
				solutionAtDataTimes[count] = firstSolver.getSolution();
				count++;
			}	
		}
		oldLike += evalLogLikelihood(data, solutionAtDataTimes, size, varData);
	}
	oldLike /= internalMC; 
	double oldPrior = evalLogPrior(paramVec, priorMean, priorVariance, nParam);
	
	VectorXd paramVecN(nParam); 	
	std::vector<VectorXd> mcmcPath(nStepsMC, VectorXd(nParam));
	std::vector<VectorXd> solutionAtDataTimes(nData, VectorXd(size));
	VectorXd w(nParam);
	std::vector<double> paramStd(nParam);

	for (int i = 0; i < nStepsMC; i++)	
	{
		if (i % 50 == 0) {
			std::cout << "iteration: " << i << std::endl;
			std::cout << "likelihood: " << oldLike << std::endl;
			std::cout << "step: " << gamma << std::endl;
			std::cout << "prior: " << oldPrior << std::endl;
			std::cout << "param: " << paramVec.transpose() << std::endl;
		}

		for (int j = 0; j < nParam; j++) {
			w(j) = gamma * normal(wGenerator);	
		}
	
		paramVecN = paramVec + w;
		for (int j = 0; j < nParam; j++) {
			paramStd[j] = paramVecN(j);
		}
	
		double lVec[internalMC];	
		int MCindex; 
		#pragma omp parallel for num_threads(12) private(MCindex)
		for (MCindex = 0; MCindex < internalMC; MCindex++) {	
			std::vector<VectorXd> solutionAtDataTimesPar(nData, VectorXd(size));
			sProbMethod<detSROCK> solver(size, h, initialCond, paramStd, odeFunc, sigma, nLevels, damping); 		
			double t = 0; 
			int count = 0;
			while (t < finalTime) {
				solver.oneStep(solGenerator, h);
				t = t + h;
				if (std::abs(t - dataTimes[count]) < h / 10.0) {
					solutionAtDataTimesPar[count] = solver.getSolution();
					count++;
				}	
			}
			lVec[MCindex] = evalLogLikelihood(data, solutionAtDataTimesPar, size, varData);
		}
		l = 0.0;		
		for (int j = 0; j < internalMC; j++) {
			l += lVec[j];
		}
		l /= internalMC;

		// Add computations with old value of parameters
		for (int j = 0; j < nParam; j++) {
			paramStd[j] = paramVec(j);
		}

		double oldParamL = 0.0;
		#pragma omp parallel for num_threads(12) private(MCindex)
		for (MCindex = 0; MCindex < internalMC; MCindex++) {	
			std::vector<VectorXd> solutionAtDataTimesPar(nData, VectorXd(size));
			sProbMethod<detSROCK> solverOld(size, h, initialCond, paramStd, odeFunc, sigma, nLevels, damping); 		
			double t = 0; 
			int count = 0;
			while (t < finalTime) {
				solverOld.oneStep(solGenerator, h);
				t = t + h;
				if (std::abs(t - dataTimes[count]) < h / 10.0) {
					solutionAtDataTimesPar[count] = solverOld.getSolution();
					count++;
				}	
			}
			lVec[MCindex] = evalLogLikelihood(data, solutionAtDataTimesPar, size, varData);
		}

		for (int j = 0; j < internalMC; j++) {
			oldParamL += lVec[j];
		}
		oldParamL /= internalMC;
 
		prior = evalLogPrior(paramVecN, priorMean, priorVariance, nParam);
		double u = log(unif(wGenerator));
		double alpha = std::min<double>(0.0, (l + prior) - (oldParamL + oldPrior));	
//		double alpha = std::min<double>(0.0, (l + prior) - (oldLike + oldPrior));
		if (u < alpha) {
			paramVec = paramVecN;
			oldLike = l;
			oldPrior = prior;
//			gamma *= 2.0; 
			gamma = std::min<double>(gamma * 2.0, 1.0);
		} else {
//			gamma /= 2.0; 
			gamma = std::max<double>(gamma / 2.0, 0.01);
		}	
		mcmcPath[i] = paramVec;	
	}
	return mcmcPath;
}

