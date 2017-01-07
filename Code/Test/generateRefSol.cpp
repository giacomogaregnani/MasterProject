#include <cmath>
#include <vector>
#include <random>
#include <Solver.hpp>
#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include "problems.hpp"
#include <iomanip>
#include <unsupported/Eigen/MatrixFunctions>

using namespace Eigen;

#define PI 3.1415926535897

int main(int argc, char* argv[])
{
    // =========================
    // INITIALIZATION
    // =========================
    problems problem = BRUSS;

    // Set problem
    odeDef testODE;
    testODE.ode = problem;
    setProblem(&testODE);
    //

    std::vector<double> paramList = testODE.refParam;
    std::string filepath("refSolBruss.txt");
    // equispaced values from 0 to finalTime
    double finalTime = 10.0;
    unsigned int nData = 10;
    double noiseStdDev = 0.1;
    // =========================

    // =========================
    // SOLUTION GENERATION
    // =========================

	// REFERENCE SOLUTION
	double hRef = 0.0001;
	int nSteps = int (finalTime / hRef);

	// Prepare structures for data acquisition
    double timeSpacing = finalTime / nData;
	std::vector<double> times;
	for (unsigned int i = 1; i < nData + 1; i++) {
		times.push_back(timeSpacing * i);
	}

	// ERROR MODEL
	std::vector<VectorXd> data(nData, VectorXd(testODE.size));
	std::normal_distribution<double> disturb(0.0, noiseStdDev);
	std::default_random_engine generator{(unsigned int) time(NULL)};

    // Write all results
    std::string folder(DATA_PATH);
    std::string slash("/");
    std::ofstream fullResults;
    std::string fullResPath = "fullResultsLorenz.txt";
    std::string fullPath = folder + slash + fullResPath;

    if (argc > 1) {
        fullResults.open(fullPath, std::ofstream::out | std::ofstream::trunc);
    }

    // COMPUTE SOLUTION (EXACT FOR TEST1D AND POISSON (MATRIX EXPONENTIALS))
    if (testODE.ode == TEST1D) {
        for (size_t i = 0; i < nData; i++) {
            data[i](0) = exp(paramList[0] * times[i]) + disturb(generator);
        }
    } else if (testODE.ode == POISSON) {
        for (unsigned int i = 0; i < nData; i++) {
            data[i] = testODE.exactSol(paramList, times[i]);
            for (int j = 0; j < testODE.size; j++) {
                data[i](j) += disturb(generator);
            }
        }
    } else if (testODE.ode == VDPOL) {
        data[0](0) = -0.1863646254808130e1 + disturb(generator);
        data[0](1) = -0.1863646254808130e1 + disturb(generator);
    } else {
        RungeKutta refSolver(testODE, paramList);
        double time = 0;
        VectorXd tmpSol = testODE.initialCond;
        int count = 0;
        for (int i = 0; i < nSteps + 10; i++) {
            tmpSol = refSolver.oneStep(tmpSol, hRef);
            time = time + hRef;
            if (fmod(time, 0.01) < 1.1 * hRef) {
                fullResults << time << "\t" << tmpSol.transpose() << "\n";
            }
            if (std::abs(times[count] - time) < hRef / 2.0) {
                std::cout << time << std::endl;
                for (int j = 0; j < testODE.size; j++) {
                    data[count](j) = tmpSol(j) + disturb(generator);
                    if (testODE.ode == HIRES && data[count](j) < 0) {
                        data[count](j) = 0.0;
                    }
                }
                count++;
            }
        }
    }

    std::ofstream results;
    std::string finalpath = folder + slash + filepath;
    results.open(finalpath , std::ofstream::out | std::ofstream::trunc);

    results << finalTime << "\n";
	results << nData << "\n";
	for (auto it : times) {
		results << it << "\n";
	}
	for (auto it : data) {
		results << std::fixed << std::setprecision(30) << it.transpose() << "\n";
	}

	results.close();
    if (argc > 1) {
        fullResults.close();
    }
	return 0;
}