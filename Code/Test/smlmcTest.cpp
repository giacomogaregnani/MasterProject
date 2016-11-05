#include <iostream>
#include <Solver.hpp>
#include <Eigen/Dense>
#include <fstream>
#include <string>
#include "problems.hpp"

using namespace Eigen;

int main(int argc, char* argv[])
{
	std::ofstream results; 
	std::string folder(DATA_PATH);
	std::string filepath(argv[1]);
	std::string slash("/");
	std::string finalpath = folder + slash + filepath;

	results.open(finalpath , std::ofstream::out | std::ofstream::trunc);
	std::cout << finalpath << std::endl;

	std::default_random_engine generator{(unsigned int) time(NULL)};

	problems problem = TEST1D;
	VectorXd (*odeFunc) (VectorXd, std::vector<double>&);
	int size;
	VectorXd initialCond;
	setProblem(&initialCond, &odeFunc, problem, &size);

	double hRef = 0.00001;
	double finalTime = 20.0;
	int nSteps = int (finalTime / hRef);
	
	// LORENZ
	//std::vector<double> parameters = {10.0, 28.0, 8.0 / 3.0};
	// TEST1D
	std::vector<double> parameters = {-10.0};	


	ProbMethod<EulerForward> refSolver(size, hRef, initialCond, parameters, odeFunc, 0.0);
	for (int i = 0; i < nSteps; i++){
		refSolver.oneStep(generator, hRef);
	}
	VectorXd solution(size);
	solution = refSolver.getSolution(); 

	// s-MLMC
	std::vector<double> desiredAccuracy = {0.1, 0.01, 0.001};//, 0.0001};
	std::vector<long int> cost;	
	VectorXd result;

	int nExp = 1;
	double finalErr;		
	 
	double rho = 10;
	double damping = 0.1;
	
	for (auto it : desiredAccuracy) {
		finalErr = 0.0;
		for (int nMCMC = 0; nMCMC < nExp; nMCMC++) {
			sMLMC<detSROCK> multilevel(size, it, initialCond, parameters, odeFunc, 0.05, 1, finalTime, damping, rho);
			if (nMCMC == 0) cost.push_back(multilevel.cost());
			std::cout << "==================" << std::endl;
			std::cout << "Computational cost: " << cost.back() << std::endl;
			result = multilevel.compute();
			std::cout << "Desired accuracy: " << it << std::endl;
			std::cout << "==================" << std::endl;
			std::cout << "Reference Solution: " << std::endl;
			std::cout << solution << std::endl;
			std::cout << "MLMC Estimation: " << std::endl;	
			std::cout << result << std::endl;
			double error = std::abs(solution.norm() - result.norm());
			std::cout << "Error = " << error << std::endl; 
			finalErr += error;
		}
		results << finalErr / nExp << "\t" << cost.back() << "\n";
	}

	// MLMC
	std::vector<long int> costNonS;	

	for (auto it : desiredAccuracy) {
		finalErr = 0.0;
		for (int nMCMC = 0; nMCMC < nExp; nMCMC++) {
			MLMC<EulerForward> multilevel(size, it, initialCond, parameters, odeFunc, 0.05, 1, finalTime);
			if (nMCMC == 0) costNonS.push_back(multilevel.cost());
			std::cout << "==================" << std::endl;
			std::cout << "Computational cost: " << costNonS.back() << std::endl;
			result = multilevel.compute();
			std::cout << "Desired accuracy: " << it << std::endl;
			std::cout << "==================" << std::endl;
			std::cout << "Reference Solution: " << std::endl;
			std::cout << solution << std::endl;
			std::cout << "MLMC Estimation: " << std::endl;	
			std::cout << result << std::endl;
			double error = (solution - result).norm();
			std::cout << "Error = " << error << std::endl; 
			finalErr += error;
		}
		results << finalErr / nExp << "\t" << costNonS.back() << "\n";
	}
	// CLOSE FILE
	results.close();
	return 0;
}

