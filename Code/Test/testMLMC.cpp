#include <cmath>
#include <vector>
#include <random>
#include <Solver.hpp>
#include <Eigen/Dense>
#include <ctime>
#include <iostream>
#include <fstream>
#include "problems.hpp"

using namespace Eigen;

#define PI 3.1415926535897
#define M_EPS 0.0000001

double phi(VectorXd, VectorXd, double);

int main(int argc, char* argv[])
{
   // Choose parameter list (on which inference will be done)
    std::vector<double> paramList = {0.2, 0.2, 3.0};

    // PROBLEM DATA
    problems problem = FITZNAG;
    VectorXd (*odeFunc) (VectorXd, std::vector<double>&);
    int size;
    VectorXd initialCond;
    setProblem(&initialCond, &odeFunc, problem, &size);

    // DATA ACQUISITION FROM REFSOL
    std::string refFile("refSolFitznag.txt");
    std::fstream refSolution;
    std::string slash("/");
    std::string folder(DATA_PATH);
    refSolution.open(folder + slash + refFile, std::ios_base::in);
    double finalTime = 0.4;
    VectorXd reference(size);
    for (int j = 0; j < size; j++) {
        refSolution >> reference(j);
    }

    // SET THE SMALLEST TIMESTEP, INITIALIZE A MLMC AND SOLVE
    double hL = 0.01;
    double sigma = 0.5;
    MLMC<EulerForward> solver(size, 0.1, initialCond, paramList,
                              odeFunc, sigma, 1, finalTime, true, hL,
                              &phi, VectorXd(), 0.0);
    double result = solver.compute();

    // REFERENCE SOLUTION AND ERROR COMPUTATION
    double refResult = reference.dot(reference);
    double error = std::abs(result - refResult);
    std::cout << "reference solution: " << refResult << std::endl
              << "MLMC estimation: " << result << std::endl
              << "error: " << error << std::endl;

    return 0;
}

double phi(VectorXd unused, VectorXd x, double unusedVar)
{
    return x.dot(x);
}
