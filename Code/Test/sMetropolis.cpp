#include <cmath>
#include <random>
#include <iostream>
#include "mcmcTools.hpp"

#define PI 3.1415926535897	

std::vector<VectorXd> sMetropolisHastings(odeDef odeModel, std::vector<double>& param, double sigma,
										  double h, double finalTime, std::vector<VectorXd>& data, std::vector<double>& dataTimes,
										  VectorXd priorMean, VectorXd priorVariance, int internalMC,
										  int nStepsMC, double damping, double varData)
{
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

	double stiffIndex;
	stiffIndex = std::abs(powerMethod(odeModel.odeJac(odeModel.initialCond, param), 1.0));
	int nLevels = static_cast<int>(std::max(2.0, std::ceil(sqrt(0.5 * stiffIndex * h))));

	// Evaluate posterior for initial guess
	sProbMethod<RKC> firstSolver(odeModel.size, h, odeModel.initialCond, param, odeModel.odeFunc, sigma, nLevels, damping); 		
	for (int j = 0; j < internalMC; j++) {	
		double t = 0; 
		int count = 0;
		std::vector<VectorXd> solutionAtDataTimes(nData, VectorXd(odeModel.size));
		while (t < finalTime) {
			firstSolver.oneStep(solGenerator, h);
			t = t + h;
			if (std::abs(t - dataTimes[count]) < h / 10.0) {
				std::cout << t << std::endl;
				solutionAtDataTimes[count] = firstSolver.getSolution();
				std::cout << "solution: " << solutionAtDataTimes[count].transpose() << std::endl;
				count++;
			}
		}
		oldLike += evalLogLikelihood(data, solutionAtDataTimes, odeModel.size, varData);
	}
	oldLike /= internalMC; 
	double oldPrior = evalLogPrior(paramVec, priorMean, priorVariance, nParam);
	
	VectorXd paramVecN(nParam); 	
	std::vector<VectorXd> mcmcPath(nStepsMC, VectorXd(nParam));
	std::vector<VectorXd> solutionAtDataTimes(nData, VectorXd(odeModel.size));
	VectorXd w(nParam);
	std::vector<double> paramStd(nParam);


	// Initialize RAM
	double gamma = 1e-4;
	double desiredAlpha = 0.25;
	MatrixXd S = RAMinit(gamma, desiredAlpha, nParam);

	for (int i = 0; i < nStepsMC; i++)	
	{
		if (i % 50 == 0) {
			std::cout << "iteration: " << i << std::endl;
			std::cout << "likelihood: " << oldLike << std::endl;
			std::cout << "step: " << S << std::endl;
			std::cout << "prior: " << oldPrior << std::endl;
			std::cout << "param: " << paramVec.transpose() << std::endl;
		}

		for (int j = 0; j < nParam; j++) {
			w(j) = normal(wGenerator);
		}
		paramVecN = paramVec + S * w;
		for (int j = 0; j < nParam; j++) {
			paramStd[j] = paramVecN(j);
		}

		stiffIndex = 1.2 * std::abs(powerMethod(odeModel.odeJac(odeModel.initialCond, paramStd), 1.0));
		nLevels = static_cast<int>(std::max(2.0, std::ceil(sqrt(0.5 * stiffIndex * h))));

		double lVec[internalMC];	
		int MCindex; 
		#pragma omp parallel for num_threads(24) private(MCindex)
		for (MCindex = 0; MCindex < internalMC; MCindex++) {	
			std::vector<VectorXd> solutionAtDataTimesPar(nData, VectorXd(odeModel.size));
			sProbMethod<RKC> solver(odeModel.size, h, odeModel.initialCond, paramStd, odeModel.odeFunc, sigma, nLevels, damping); 		
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
			lVec[MCindex] = evalLogLikelihood(data, solutionAtDataTimesPar, odeModel.size, varData);
		}
		l = 0.0;		
		for (int j = 0; j < internalMC; j++) {
			l += lVec[j];
		}
		l /= internalMC;
		std::cout << "param = " << paramVecN.transpose()
				  << " stiffness = " << stiffIndex
				  << " nLevels = " << nLevels
				  << " likelihood = " << l
				  << std::endl;

		// Add computations with old value of parameters
		for (int j = 0; j < nParam; j++) {
			paramStd[j] = paramVec(j);
		}

		stiffIndex = 1.2 * std::abs(powerMethod(odeModel.odeJac(odeModel.initialCond, paramStd), 1.0));
		nLevels = static_cast<int>(std::max(2.0, std::ceil(sqrt(0.5 * stiffIndex * h))));

		double oldParamL = 0.0;
		#pragma omp parallel for num_threads(24) private(MCindex)
		for (MCindex = 0; MCindex < internalMC; MCindex++) {	
			std::vector<VectorXd> solutionAtDataTimesPar(nData, VectorXd(odeModel.size));
			sProbMethod<RKC> solverOld(odeModel.size, h, odeModel.initialCond, paramStd, odeModel.odeFunc, sigma, nLevels, damping); 		
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
			lVec[MCindex] = evalLogLikelihood(data, solutionAtDataTimesPar, odeModel.size, varData);
		}

		for (int j = 0; j < internalMC; j++) {
			oldParamL += lVec[j];
		}
		oldParamL /= internalMC;
 
		prior = evalLogPrior(paramVecN, priorMean, priorVariance, nParam);
		double u = log(unif(wGenerator));
		double alpha = std::min<double>(0.0, (l + prior) - (oldParamL + oldPrior));	
		if (u < alpha) {
			paramVec = paramVecN;
			oldLike = l;
			oldPrior = prior;
		} 
		mcmcPath[i] = paramVec;

		// Update RAM
		S = RAMupdate(S, w, exp(alpha), desiredAlpha, nParam, i + 1);
	}
	return mcmcPath;
}

