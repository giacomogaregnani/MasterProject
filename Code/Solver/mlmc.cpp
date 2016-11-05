#include "Solver.hpp"

// MULTI LEVEL
template <class T>
MLMC<T>::MLMC(int n, double accuracy,
      VectorXd initialCondition, std::vector<double> parameters,
      VectorXd (*func) (VectorXd, std::vector<double>&),
      double stoch, int order, double finalTime, bool fromTimestep, double hL,
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
		nLevels = 1 + static_cast<int>(ceil(1.0 / static_cast<double>(order) * std::abs(log2(accuracy))));
	} else {
		nLevels = 1 + static_cast<int>(ceil(log2(finalTime / hL)));
	}

	// Computation of timesteps, number of steps and number of trajectories per level	
	for (int i = nLevels - 1; i > -1; i--){	
		nSteps.push_back(static_cast<int>(pow(2.0, i)));
		timesteps.push_back(finalTime / nSteps.back());
	}
	for (auto it : timesteps) {
		nTrajectories.push_back(static_cast<long int>(pow(it, 2 * order) * static_cast<double>(nLevels) /
                                                      pow(timesteps.front(), 2 * order)));
		solver.push_back(std::make_shared<ProbMethod<T>>(n, it, initialCond, parameters, func, stoch));	
	}

	phiMLMC = phi;
	dataMCMC = data;
	varDataMCMC = varianceData;

	/* std::cout << "======================" << std::endl;
	std::cout << "Initialized a MLMC estimator" << std::endl
		  << "Number of levels: " << nLevels << std::endl;

	for (int i = 0; i < nLevels; i++){
		std::cout << "======================" << std::endl;
		std::cout << "Level " << i << std::endl
			  << "timestep: " << timesteps[i] << std::endl
			  << "N. of trajectories: " << nTrajectories[i] << std::endl;
	} */
}

template <class T>
long int MLMC<T>::cost(void)
{
	long int cost = 0;
	for (int i = 0; i < nLevels; i++) {
		cost += nTrajectories[i] * nSteps[i];
	}	
	return cost;
}

template <class T>
double MLMC<T>::compute(void)
{
	// Initialize the level difference and the phis to zero
    // phis has one entry more as by definition phi_{-1} = 0
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
                    tmp = solver[k]->oneStepGiven(brownianPath[steps], timesteps[k]);
                }
                phis[k] = phiMLMC(dataMCMC, tmp, varDataMCMC);
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
double MLMC<T>::gethL(void)
{
	return timesteps.front();
} 

template class MLMC<RungeKutta>;
template class MLMC<EulerForward>;

