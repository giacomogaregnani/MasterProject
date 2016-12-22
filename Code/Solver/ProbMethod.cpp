#include "Solver.hpp"

template <class T> 
ProbMethod<T>::ProbMethod(odeDef ODE, double timestep, std::vector<double> paramVec,
						  double stoch)
{
	detSolver = std::make_shared<T>(ODE, paramVec);
	solution = ODE.initialCond;
	sigma = stoch;
	rootsigma = sqrt(sigma); 
	size = ODE.size;
	h = timestep;
	int detOrder = detSolver->getOrder();
	hfunc = pow(h, detOrder + 0.5);
	IC = ODE.initialCond;
	parameters = paramVec;
}

template <class T>
VectorXd& ProbMethod<T>::getSolution(void)
{
	return solution;
}

template <class T> 
void ProbMethod<T>::oneStep(std::default_random_engine& generator, double step) 
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
VectorXd ProbMethod<T>::oneStepGiven(VectorXd& noise, double step)
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
void ProbMethod<T>::resetIC(void)
{
	solution = IC;
}

template class ProbMethod<RungeKutta>;
template class ProbMethod<EulerForward>;
template class ProbMethod<EulerBackwards>;
template class ProbMethod<MidPoint>;