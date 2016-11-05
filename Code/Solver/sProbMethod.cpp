#include "Solver.hpp"

// STABILIZED PROBABILISTIC METHOD
template <class T> 
sProbMethod<T>::sProbMethod(int n, double timestep,
                          VectorXd initialCond, std::vector<double> paramVec,
                          VectorXd (*func) (VectorXd, std::vector<double>&),
	                  double stoch, int nRKCStages, double damping)
{
	detSolver = std::make_shared<T>(n, func, paramVec, nRKCStages, damping);
	detSolver->computeStageCoeff(timestep);
	solution = initialCond;
	sigma = stoch;
	rootsigma = sqrt(sigma); 
	size = n;
	h = timestep;
	int detOrder = detSolver->getOrder();
	hfunc = pow(h, detOrder + 0.5);
	IC = initialCond;
	parameters = paramVec;
	nStages = nRKCStages;
}

template <class T>
VectorXd& sProbMethod<T>::getSolution(void)
{
	return solution;
}

template <class T> 
void sProbMethod<T>::oneStep(std::default_random_engine& generator, double step) 
{
	VectorXd noise(size);
	for (int i = 0; i < size; i++){	
		noise(i) = normalDist(generator);
	}

	if (step == h) {
		VectorXd detSol = detSolver->oneStep(solution, h);
		solution = detSol + hfunc * rootsigma * noise;
	} else {
		double hfunctmp = pow(h, detSolver->getOrder() + 0.5);
		VectorXd detSol = detSolver->oneStep(solution, step);
		solution = detSol + hfunctmp * rootsigma * noise;
	}
}

template <class T>
VectorXd sProbMethod<T>::oneStepGiven(VectorXd& noise, double step)
{
	if (step == h) {
		VectorXd detSol = detSolver->oneStep(solution, h);
		solution = detSol + hfunc * rootsigma * noise;
	} else {
		double hfunctmp = pow(h, detSolver->getOrder() + 0.5);
		VectorXd detSol = detSolver->oneStep(solution, step);
		solution = detSol + hfunctmp * rootsigma * noise;
	}
	return solution;
}

template <class T>
void sProbMethod<T>::resetIC(void)
{
	solution = IC;
}

template class sProbMethod<detSROCK>;
