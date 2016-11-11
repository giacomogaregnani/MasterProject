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
    problems problem = BRUSS;
    std::vector<double> paramList = {1.0};
    std::string filepath("refSolBruss.txt");
    // equispaced values from 0 to finalTime
    double finalTime = 10;
    unsigned int nData = 1;
    // =========================

    // =========================
    // KEEP THE REST
    // =========================

   	// Set problem
    odeDef testODE;
	testODE.ode = problem;
	setProblem(&testODE);

	// REFERENCE SOLUTION
	double hRef = 0.000001;
	int nSteps = int (finalTime / hRef);

	// Prepare structures for data acquisition
    double timeSpacing = finalTime / nData;
	std::vector<double> times;
	for (unsigned int i = 1; i < nData + 1; i++) {
		times.push_back(timeSpacing * i);
	}

	// ERROR MODEL
	std::vector<VectorXd> data(nData, VectorXd(testODE.size));
	std::normal_distribution<double> disturb(0.0, 1e-1);
	std::default_random_engine generator{(unsigned int) time(NULL)};

    // COMPUTE SOLUTION (EXACT FOR TEST1D AND POISSON (MATRIX EXPONENTIALS))
    if (testODE.ode == TEST1D) {
        for (size_t i = 0; i < nData; i++) {
            data[i](0) = exp(paramList[0] * times[i]) + disturb(generator);
        }
    } else if (testODE.ode == POISSON) {
        MatrixXd A(testODE.size, testODE.size);
        for (int i = 0; i < testODE.size - 1; i++) {
            A(i, i) = -2.0;
            A(i, i + 1) = 1.0;
            A(i + 1, i) = 1.0;
        }
        A(testODE.size - 1, testODE.size - 1) = -2.0;

        for (unsigned int i = 0; i < nData; i++) {
            MatrixXd B = paramList[0] * times[i] * A;
            data[i] = B.exp() * testODE.initialCond;
            for (int j = 0; j < testODE.size; j++) {
                data[i](j) += disturb(generator);
            }
        }
    } else {
        ProbMethod<EulerForward> refSolver(testODE.size, hRef, testODE.initialCond, paramList, testODE.odeFunc, 0.0);
        double time = 0;
        VectorXd tmpSol(testODE.size);
        int count = 0;
        for (int i = 0; i < nSteps; i++) {
            refSolver.oneStep(generator, hRef);
            time = time + hRef;
            if (std::abs(times[count] - time) < hRef / 10.0) {
                std::cout << time << std::endl;
                tmpSol = refSolver.getSolution();
                for (int j = 0; j < testODE.size; j++) {
                    data[count](j) = tmpSol(j) + disturb(generator);
                }
                count++;
            }
        }
    }

    std::ofstream results;
    std::string folder(DATA_PATH);
    std::string slash("/");
    std::string finalpath = folder + slash + filepath;

    results.open(finalpath , std::ofstream::out | std::ofstream::trunc);
    std::cout << finalpath << std::endl;

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