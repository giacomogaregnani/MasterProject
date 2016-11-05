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
	
	double sigma = 0.1;

	VectorXd (*odeFunc) (VectorXd, std::vector<double>&);
	odeFunc = &fitznag;
	std::vector<double> parameters = {0.2, 0.2, 3.0};

	int nStages = 5;
	double damping = 5.9;
	
	int size = 2;
	
	double hRef = 0.1;
	double finalTime = 20.0;
	int nSteps = int (finalTime / hRef);
	VectorXd solution(size);

	// generator of random numbers
	std::default_random_engine generator{(unsigned int) time(NULL)};
	
	sProbMethod<detSROCK> testSolver(size, hRef, IC, parameters, odeFunc, sigma, nStages, damping);
	for (int i = 0; i < nSteps; i++){
		testSolver.oneStep(generator, hRef);
		results << testSolver.getSolution().transpose() << "\n";
	}
	solution = testSolver.getSolution();
	results << solution.transpose() << "\n"; 
	std::cout << solution.transpose() << std::endl;
	
	ProbMethod<EulerForward> refSolver(size, hRef, IC, parameters, odeFunc, sigma);
	for (int i = 0; i < nSteps; i++){
		refSolver.oneStep(generator, hRef);
		results << refSolver.getSolution().transpose() << "\n";	
	}	
	solution = refSolver.getSolution(); 
	results << solution.transpose() << "\n";
	std::cout << solution.transpose() << std::endl;

	results.close();
	return 0;
}
