#include <iostream>
#include <Solver.hpp>
#include <Eigen/Dense>
#include <fstream>
#include <string>
#include <cmath>
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
	VectorXd (*odeFunc) (VectorXd, std::vector<double>&);
	int size;
	VectorXd initialCond;
	problems problem = BRUSS;
	setProblem(&initialCond, &odeFunc, problem, &size);

	double hRef = 0.00001;
	double alpha = 1.0 / 50.0;
	double dX = 1.0 / (static_cast<double>(size / 2.0 + 1.0));
	std::vector<double> parameters = {alpha, dX};

	// interesting times definition
	std::vector<double> times = {1.0, 2.0, 3.0, 4.0, 10.0};
	double finalTime = times.back();
	int nSteps = int (finalTime / hRef);

	// Compute the solution at all interesting times	
	ProbMethod<EulerForward> refSolver(size, hRef, initialCond, parameters, odeFunc, 0.0);
	//sProbMethod<detSROCK> refSolver(size, hRef, initialCond, parameters, odeFunc, 0.0, 200, 5.0);
	double tRef = 0;
	VectorXd refX(size);	
	std::vector<double> refSolution = {};
	int count = 0;	
	for (int i = 0; i < nSteps + 1; i++){
		refSolver.oneStep(generator, hRef);
		tRef += hRef;
		if (std::abs(tRef - times[count]) < hRef / 10.0) {
			refX = refSolver.getSolution();
			refSolution.push_back(refX.dot(refX)); 
			count++;
			std::cout << tRef << std::endl;
		}	
	}

	// MLMC
	std::vector<double> desiredAccuracy = {0.01};//, 0.01, 0.001, 0.0001};
	std::vector<long int> cost;	
	double MLMCresult;	

	int nExp = 1;
	double finalErr;		
	double hL;
//	MLMC<EulerForward>* multilevel = nullptr;
	sMLMC<detSROCK>* multilevel = nullptr;
	double damping = 0.1;
	double rho = 300.0;
	
	for (auto it : desiredAccuracy) {
		for (size_t nTime = 0; nTime < times.size(); nTime++) {
			finalErr = 0.0;
			for (int nMCMC = 0; nMCMC < nExp; nMCMC++) {
				if (nTime == 0) {
					//multilevel = new MLMC<EulerForward>(size, it, initialCond, parameters, 
					//		      odeFunc, 0.1, 1, times[0], false, 0.0);
					multilevel = new sMLMC<detSROCK>(size, it, initialCond, parameters,
							      odeFunc, 0.1, 1, times[0], damping, rho, false, 0.0); 
					hL = multilevel->gethL();	
				} else {
					//multilevel = new MLMC<EulerForward>(size, it, initialCond, parameters, 
					//		      odeFunc, 0.1, 1, times[nTime], true, hL);
					multilevel = new sMLMC<detSROCK>(size, it, initialCond, parameters,
							      odeFunc, 0.1, 1, times[nTime], damping, rho, true, hL); 
				}
				if (nMCMC == 0) cost.push_back(multilevel->cost());
				std::cout << "==================" << std::endl;
				std::cout << "Computational cost: " << cost.back() << std::endl;
				MLMCresult = multilevel->compute();
				std::cout << "Desired accuracy: " << it << std::endl;
				std::cout << "==================" << std::endl;
				std::cout << "Reference Solution: " << std::endl;
				std::cout << refSolution[nTime] << std::endl;
				std::cout << "MLMC Estimation: " << std::endl;	
				std::cout << MLMCresult << std::endl;
				double error = std::abs(refSolution[nTime] - MLMCresult);
				std::cout << "Error = " << error << std::endl; 
				finalErr += error;
			}
			results << finalErr / nExp << "\t" << cost.back() << "\n";
			delete multilevel;
		}	
	}
	
	// CLOSE FILE
	results.close();
	return 0;
}

