#include <Solver.hpp>
#include <iostream>
#include <Eigen/Dense>
#include <fstream>
#include <vector>
#include <memory>
#include <cmath>
#include <random>
#include <ctime>

using namespace Eigen;

VectorXd lorenz(VectorXd argument);
VectorXd test(VectorXd argument);
VectorXd fitznag(VectorXd argument);

enum problems {
	FITZNAG,
	LORENZ,
	TEST
};

int main(void)
{
	std::ofstream results; 
	results.open(DATA_PATH "/orders.txt", std::ofstream::out | std::ofstream::trunc);

	VectorXd (*odeFunc) (VectorXd);
	problems problem = FITZNAG;

	int size;
	VectorXd initialCond;
	std::default_random_engine generator{(unsigned int) time(NULL)};

	switch (problem) {
		case FITZNAG:
			size = 2;
			initialCond.resize(size);
			initialCond(0) = -1.0;
			initialCond(1) = 1.0;
			odeFunc = &fitznag;
			break;
		case LORENZ:
			size = 3;
			initialCond.resize(size);
			initialCond(0) = -10.0;
			initialCond(1) = -1.0;
			initialCond(2) = 40.0;
			odeFunc = &lorenz;
			break;	
		case TEST:
			size = 1;
			initialCond.resize(size);
			initialCond(0) = 1.0;
			odeFunc = &test;
			break;
	}

	// Reference solution

	double hRef = 0.00001;
	double finalTime = 10.0;
	int nSteps = int (finalTime / hRef);
	
	ProbMethod<RungeKutta> refSolver(size, hRef, initialCond, odeFunc, 0.0);
	for (int i = 0; i < nSteps; i++){
		refSolver.oneStep(generator, hRef);
	}
	VectorXd solution(size);
	solution = refSolver.getSolution(); 
	results << solution << "\n";

	// MonteCarlo

	int nMC = 100000;
	std::vector<double> h(4);
	h[0] = 0.25;
	for (unsigned int i = 1; i < h.size(); i++) {
		h[i] = h[i - 1] / 2.0; 
	}
	double stoch = 0.1;
		
	for (auto it : h) {
		std::cout << it << std::endl;	
		for (int j = 0; j < size; j++) {
			solution(j) = 0.0;
		}
		for (int j = 0; j < nMC; j++){
			ProbMethod<RungeKutta> solver(size, it, initialCond, odeFunc, stoch);
			double t = 0;
			while (t < finalTime){
				solver.oneStep(generator, it);
				t = t + it;
			}
			solution += solver.getSolution();
		}
		solution /= nMC;
		results << solution << "\n";
	}
	
	results.close();
	return 0; 
}

VectorXd lorenz(VectorXd argument) {
	double sigma = 10, rho = 28, beta = 8.0/3.0;

	VectorXd result(3);
	result(0) = sigma * (argument(1) - argument(0));
	result(1) = argument(0) * (rho - argument(2)) - argument(1);
	result(2) = argument(0) * argument(1) - beta * argument(2);
	return result;
}

VectorXd test(VectorXd argument) {
	MatrixXd A(3, 3);
	
	for (int i = 0; i < 3; i++) {
		A(i, i) = -2;
	}
	
	return A * argument;
}

VectorXd fitznag(VectorXd argument) {
	double a = 0.2, b = 0.2, c = 3;
	
	VectorXd result(2);
	result(0) = c * (argument(0) - argument(0) * argument(0) * argument(0) / 3.0 + argument(1));
	result(1) = - 1.0 / c * (argument(0) - a + b * argument(1));

	return result;
}

