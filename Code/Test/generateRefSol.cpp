#include <cmath>
#include <vector>
#include <random>
#include <Solver.hpp>
#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include "problems.hpp"
#include <unsupported/Eigen/MatrixFunctions>

using namespace Eigen;

#define PI 3.1415926535897	

int main(int argc, char* argv[])
{
    // =========================
    // CHANGE THOSE VALUES
    // =========================
    problems problem = LORENZ;
    std::vector<double> paramList = {10.0, 28.0, 8.0 / 3.0};
    std::string filepath("refSolLorenz.txt");
    double finalTime = 20;
    unsigned int nData = 20;
    // =========================

    // =========================
    // KEEP THE REST
    // =========================
    std::ofstream results;
	std::string folder(DATA_PATH);
	std::string slash("/");
	std::string finalpath = folder + slash + filepath;

	results.open(finalpath , std::ofstream::out | std::ofstream::trunc);
	std::cout << finalpath << std::endl;
	
	// Set problem
	VectorXd (*odeFunc) (VectorXd, std::vector<double>&);
	int size;
	VectorXd initialCond;
	setProblem(&initialCond, &odeFunc, problem, &size);

	// REFERENCE SOLUTION
	double hRef = 0.000001;
	int nSteps = int (finalTime / hRef);

	// Prepare structures for data acquisition
	std::vector<double> times;
	for (unsigned int i = 1; i < nData + 1; i++) {
		times.push_back(static_cast<double>(i));
	}

	// ERROR MODEL
	std::vector<VectorXd> data(nData, VectorXd(size));	
	std::normal_distribution<double> disturb(0.0, 1e-3);
	std::default_random_engine generator{(unsigned int) time(NULL)};
	
    // COMPUTE SOLUTION (EXACT FOR TEST1D AND POISSON (MATRIX EXPONENTIALS))
    if (problem == TEST1D) {
        for (size_t i = 0; i < nData; i++) {
            data[i](0) = exp(paramList[0] * times[i]) + disturb(generator);
        }
    } else if (problem == POISSON) {
        MatrixXd A(size, size);
        for (int i = 0; i < size - 1; i++) {
            A(i, i) = -2.0;
            A(i, i + 1) = 1.0;
            A(i + 1, i) = 1.0;
        }
        A(size - 1, size - 1) = -2.0;

        for (unsigned int i = 0; i < nData; i++) {
            MatrixXd B = paramList[0] * times[i] * A;
            data[i] = B.exp() * initialCond;
            for (int j = 0; j < size; j++) {
                data[i](j) += disturb(generator);
            }
        }
    } else {
        ProbMethod<EulerForward> refSolver(size, hRef, initialCond, paramList, odeFunc, 0.0);
        double time = 0;
        VectorXd tmpSol(size);
        int count = 0;
        for (int i = 0; i < nSteps; i++) {
            refSolver.oneStep(generator, hRef);
            time = time + hRef;
            if (std::abs(times[count] - time) < hRef / 10.0) {
                std::cout << time << std::endl;
                tmpSol = refSolver.getSolution();
                for (int j = 0; j < size; j++) {
                    data[count](j) = tmpSol(j) + disturb(generator);
                }
                count++;
            }
        }
    }

	results << finalTime << "\n";
	results << nData << "\n";
	for (auto it : times) {
		results << it << "\n";
	}
	for (auto it : data) {
		results << it.transpose() << "\n";
	}
	
	results.close();
	return 0;
}