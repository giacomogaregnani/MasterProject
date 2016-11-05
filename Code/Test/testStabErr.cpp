#include <Solver.hpp>
#include <Eigen/Dense>
#include "problems.hpp"
#include <random>
#include <ctime>
#include <iostream>
#include <fstream>

using namespace Eigen;

int main(int argc, char* argv[])
{
	if (argc == 1) {
		std::cout << "give a name for result file (written in data/namefile.txt)" << std::endl;
		return 1;
	}

	std::ofstream results; 
	std::string folder(DATA_PATH);
	std::string filepath(argv[1]);
	std::string slash("/");
	std::string finalpath = folder + slash + filepath;

	results.open(finalpath , std::ofstream::out | std::ofstream::trunc);
	std::cout << finalpath << std::endl;

	VectorXd IC(2);
	IC(0) = -1.0;
	IC(1) = 1.0;
	
	double sigma = 0.0;

	VectorXd (*odeFunc) (VectorXd, std::vector<double>&);
	odeFunc = &fitznag;
	std::vector<double> parameters = {0.2, 0.2, 3.0};

	int nStages = 5;
	double damping = 5.9;
	
	int size = 2;
	
	double hRef = 0.000001;
	double finalTime = 10.0;
	int nSteps = int (finalTime / hRef);
	VectorXd solution(size);

	// generator of random numbers
	std::default_random_engine generator{(unsigned int) time(NULL)};
	
	ProbMethod<EulerForward> refSolver(size, hRef, IC, parameters, odeFunc, sigma);
	for (int i = 0; i < nSteps; i++){
		refSolver.oneStep(generator, hRef);
	}	
	solution = refSolver.getSolution(); 

	std::vector<double> h = {1, 0.5, 0.1, 0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001};
	
	std::vector<int> nStepsStab;

	for (auto it : h) {
		nStepsStab.push_back(static_cast<int> (finalTime / it));	
	}
	
	std::vector<double> err;
	VectorXd stabSol(size);

	int count = 0;
	for (auto it : h) {	
		std::cout << it << std::endl;
		sProbMethod<detSROCK> testSolver(size, it, IC, parameters, odeFunc, sigma, nStages, damping);
		for (int i = 0; i < nStepsStab[count]; i++){
			testSolver.oneStep(generator, hRef);
		}
		count++;
		stabSol = testSolver.getSolution();	
		err.push_back((stabSol - solution).norm());
	}
	
	for (auto it : err) {
		results << it << "\n";
	}
	
	results.close();
	return 0;
}
