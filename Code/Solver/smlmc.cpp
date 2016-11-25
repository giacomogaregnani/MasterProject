#include "Solver.hpp"

// MULTI LEVEL
template <class T>
sMLMC<T>::sMLMC(int n, double accuracy,
	        VectorXd initialCondition, std::vector<double> parameters,
	        VectorXd (*func) (VectorXd, std::vector<double>&),
	        double stoch, int order, double finalTime, double damping, double rho,
	        bool fromTimestep, double hL,
	        double(*phi) (VectorXd, VectorXd, double),
	        VectorXd data, double varianceData)
{		
	size = n;
	epsilon = accuracy;
	initialCond = initialCondition;

	// If the user passed a desired smallest timestep, compute the number
	// of levels accordingly. Otherwise, decide the number of levels based
	// on the required accuracy
	if (!fromTimestep) {
		nLevels = 1 + static_cast<int>(std::ceil(1.0 / static_cast<double>(order) * std::abs(log2(accuracy))));
	} else {
		nLevels = 1 + static_cast<int>(std::round(log2(finalTime / hL)));
	}

	for (int i = nLevels - 1; i > -1; i--){
		nSteps.push_back(static_cast<int>(pow(2, i)));
		timesteps.push_back(finalTime / nSteps.back());
		nStages.push_back(static_cast<int>(std::ceil(std::max(2.0, sqrt(3.0 * finalTime * rho)
																   * pow(2.0, -static_cast<double>(i) / 2)))));
	}

	int count = 0;
	for (auto it : timesteps) {	
		nTrajectories.push_back(static_cast<long int>(std::ceil(pow(it, 2 * order) * static_cast<double>(nLevels)
																/ pow(timesteps.front(), 2 * order))));
		solver.push_back(std::make_shared<sProbMethod<T>>(n, it, initialCond, parameters, func, stoch, nStages[count++], damping));
	}

	phiMLMC = phi;
	dataMC = data;
	varData = varianceData;
}

template <class T>
double sMLMC<T>::compute(void)
{
	// Initialize the level difference and the phis to zero
	std::vector<double> phis(nLevels + 1, 0.0);
	std::vector<double> levelDiff(nLevels, 0.0);
	VectorXd tmp(size);
	int indexJump;
	
	for (int i = 0; i < nLevels; i++) {
		std::vector<VectorXd> brownianPath(nSteps[i], VectorXd(size));
		long int simTrajectories;
		if (i == 0) {
			simTrajectories = nTrajectories[i];
		} else {
			simTrajectories = nTrajectories[i] - nTrajectories[i - 1];
		}
		for (int j = 0; j < simTrajectories; j++) {
			// Create the Brownian path
			for (int k = 0; k < nSteps[i]; k++) {
				for (int sizeIt = 0; sizeIt < size; sizeIt++) {
					brownianPath[k](sizeIt) = normalDist(generator);
				}	
			}
			// Solve for all levels for which it is pertinent
			for (int k = nLevels - 1; k > i-1; k--) {
				indexJump = nSteps[i] / nSteps[k];
				for (int steps = indexJump - 1; steps < nSteps[i]; steps = steps + indexJump) {
					tmp = solver[k]->oneStepGiven(brownianPath[steps], timesteps[k], 1);
				}
				phis[k] = phiMLMC(dataMC, tmp, varData);
				solver[k]->resetIC(); 
				levelDiff[k] += phis[k] - phis[k + 1];
			}
		}
	}	

	double result = 0.0;
	for (int i = 0; i < nLevels; i++) {
		levelDiff[i] /= nTrajectories[i];		
		result += levelDiff[i];
	}

	return result;	
}

template <class T>
long int sMLMC<T>::cost(void)
{
	long int cost = 0;
	for (int i = 0; i < nLevels; i++) {
		cost += nTrajectories[i] * nSteps[i] * nStages[i];
	}	
	return cost;
}

template <class T>
double sMLMC<T>::gethL(void)
{
	return timesteps.front();
}


template class sMLMC<detSROCK>;
template class sMLMC<RKC>;
