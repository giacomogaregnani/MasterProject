#include <cmath>
#include <vector>
#include <random>
#include <Solver.hpp>
#include <fstream>
#include "problems.hpp"

using namespace Eigen;

#define PI 3.1415926535897

int main(int argc, char* argv[]) {
    // =========================
    // INITIALIZATION
    // =========================
    problems problem = HIRES;
    // Set problem
    odeDef testODE;
    testODE.ode = problem;
    setProblem(&testODE);
    //
    std::vector<double> paramList = testODE.refParam;
    std::string filepath("refSolHires.txt");
    // equispaced values from 0 to finalTime
    double finalTime = 10.0;
    // =========================

    // =========================
    // Generate a set of trajectories
    // =========================

    // REFERENCE SOLUTION
    double hRef = 0.0001;
    int nSteps = int(finalTime / hRef);

    // random generator
    std::default_random_engine generator{(unsigned int) time(NULL)};

    // Write all results
    std::string folder(DATA_PATH);
    std::string slash("/");
    std::ofstream fullResults;
    std::string fullResPath = "jacHires.txt";
    std::string fullPath = folder + slash + fullResPath;
    fullResults.open(fullPath, std::ofstream::out | std::ofstream::trunc);

    ProbMethod<EulerForward> refSolver(testODE, hRef, paramList, 0.5);
    int nRealizations = 1;

    for (int j = 0; j < nRealizations; j++) {
        for (int i = 0; i < nSteps; i++) {
            refSolver.oneStep(generator, hRef);
            if (i % 1 == 0) {
                // fullResults << refSolver.getSolution().transpose() << "\n";
                // lambda = powerMethod(, testODE.refParam), 1.0, 1000);
                fullResults << testODE.odeJac(refSolver.getSolution(), testODE.refParam)  << "\n";
            }
        }
        refSolver.resetIC();
    }

    fullResults.close();
    return 0;
}