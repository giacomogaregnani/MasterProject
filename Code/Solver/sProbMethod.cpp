#include "Solver.hpp"

// STABILIZED PROBABILISTIC METHOD
template <class T> 
sProbMethod<T>::sProbMethod(odeDef ODE, double timestep,
							std::vector<double> paramVec, double stoch,
							int nRKCStages, double damping)
{
	detSolver = std::make_shared<T>(ODE, paramVec, nRKCStages, damping);
	detSolver->computeStageCoeff(timestep);
	solution = ODE.initialCond;
	sigma = stoch;
	rootsigma = sqrt(sigma); 
	size = ODE.size;
	h = timestep;
	int detOrder = detSolver->getOrder();
	hfunc = pow(h, detOrder + 0.5);
	IC = ODE.initialCond;
	parameters = paramVec;
	nStages = nRKCStages;
}

template <class T>
VectorXd& sProbMethod<T>::getSolution(void)
{
	return solution;
}

template <class T> 
void sProbMethod<T>::oneStep(std::default_random_engine& generator, double step, int localStages)
{
	VectorXd noise(size);
	for (int i = 0; i < size; i++){	
		noise(i) = normalDist(generator);
	}

	if (step == h) {
		VectorXd detSol = detSolver->oneStep(solution, h, localStages);
		solution = detSol + hfunc * rootsigma * noise;
	} else {
		double hfunctmp = pow(h, detSolver->getOrder() + 0.5);
		VectorXd detSol = detSolver->oneStep(solution, step, localStages);
		solution = detSol + hfunctmp * rootsigma * noise;
	}
}

template <class T>
VectorXd sProbMethod<T>::oneStepGiven(VectorXd& noise, double step, int localStages)
{
	if (step == h) {
		VectorXd detSol = detSolver->oneStep(solution, h, localStages);
		solution = detSol + hfunc * rootsigma * noise;
	} else {
		double hfunctmp = pow(h, detSolver->getOrder() + 0.5);
		VectorXd detSol = detSolver->oneStep(solution, step, localStages);
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
template class sProbMethod<RKC>;