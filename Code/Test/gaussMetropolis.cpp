#include <cmath>
#include <random>
#include <iostream>
#include "mcmcTools.hpp"

#define PI 3.1415926535897	

MatrixXd diffusion(VectorXd x, std::vector<double>& param, double sigma, double h)
{
	long int size = x.size();
	MatrixXd M = MatrixXd::Identity(size, size) * h * sigma;
	return M;
}

std::vector<VectorXd> gaussMetropolisHastings(odeDef odeModel, std::vector<double>& param, double sigma,
											  double h, double finalTime, std::vector<VectorXd>& data,
											  std::vector<double>& dataTimes,
											  VectorXd priorMean, VectorXd priorVariance, int nStepsMC,
											  double varData, double* accRatio, bool isStable, StabValues stabParam)
{
    // Initialization of the acceptance ratio
	*accRatio = 0.0;

    // likelihoods and priors
	double l, prior, oldLike, oldPrior;

    // Create the random seeds
    std::default_random_engine wGenerator{(unsigned int) time(NULL)};
	std::normal_distribution<double> normal(0.0, 1.0);
	std::uniform_real_distribution<double> unif;

    //  Determine the number of data
	size_t nData = data.size();

    // Create first parameter guess
	int nParam = static_cast<int> (param.size());
	VectorXd paramVec(nParam);
	for (int i = 0; i < nParam; i++) {
		paramVec(i) = param[i];
	} 

	// compute number of step for going from data to data with time step h
	std::vector<int> nSteps(nData);
	nSteps[0] = static_cast<int>(round(dataTimes[0] / h));
	for (size_t i = 1; i < nData; i++) {
		nSteps[i] = static_cast<int>(round((dataTimes[i] - dataTimes[i - 1]) / h));
	}

    // Initialize the variance of the initial condition (zeros if the initial condition is deterministic)
    MatrixXd initialVar = varData * MatrixXd::Identity(odeModel.size, odeModel.size);

    // Compute likelihood and prior first guess of the parameters
	oldLike = 0.0;
	ThirdOrderGauss solver(odeModel, initialVar, param, &diffusion, varData, h, sigma, isStable, &stabParam);
	for (size_t i = 0; i < nData; i++) {
		oldLike += solver.oneStep(data[i], nSteps[i]);
	}
	oldPrior = evalLogPrior(paramVec, priorMean, priorVariance, nParam);

    // Initialize some support structure
	VectorXd paramVecN(nParam);
	std::vector<VectorXd> mcmcPath(nStepsMC, VectorXd(nParam));
	std::vector<double> paramStd(nParam);

    // Initialize the proposal vector
    VectorXd w(nParam);

    // Initialize RAM
    double gamma = 1e-4;
    double desiredAlpha = 0.25;
    MatrixXd S = RAMinit(gamma, desiredAlpha, nParam);

	// Only for positive parameters
	// double oldGauss, newGauss;

	for (int i = 0; i < nStepsMC; i++) {
        // Print results every 50 iterations
		if (i % 1 == 0) {
			std::cout << "iteration: " << i << std::endl;
			std::cout << "likelihood: " << oldLike << std::endl;
			std::cout << "S matrix: " << std::endl << S << std::endl;
			std::cout << "prior: " << oldPrior << std::endl;
			std::cout << "param: " << paramVec.transpose() << std::endl;
		}

		// Generate new guess for parameter
		/* bool isPositive = false;

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
		}*/

		for (int j = 0; j < nParam; j++) {
			w(j) = normal(wGenerator);
		}
		paramVecN = paramVec + S * w;

        for (int j = 0; j < nParam; j++) {
			paramStd[j] = paramVecN(j);
		}

        // Compute the likelihood and the prior for this new guess
        l = 0.0;
        ThirdOrderGauss newSolver(odeModel, initialVar, paramStd, &diffusion, varData, h, sigma, isStable, &stabParam);
		for (size_t j = 0; j < nData; j++) {
			l += newSolver.oneStep(data[j], nSteps[j]);
		}
		prior = evalLogPrior(paramVecN, priorMean, priorVariance, nParam);

		// Generate probability and update
		double u = log(unif(wGenerator));
		// double alpha = std::min<double>(0.0, l + prior + oldGauss - newGauss - oldLike - oldPrior);
        double alpha = std::min<double>(0.0, l + prior - oldLike - oldPrior);
		if (u < alpha) {
			*accRatio += 1.0 / nStepsMC;
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