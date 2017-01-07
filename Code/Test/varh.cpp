#include <cmath>
#include <vector>
#include <random>
#include <Solver.hpp>
#include <fstream>
#include "problems.hpp"
#include <iomanip>

using namespace Eigen;
using namespace std;

#define PI 3.1415926535897

int main(int argc, char* argv[]) {
    // =========================
    // INITIALIZATION
    // =========================
    problems problem = LORENZ;
    odeDef testODE;
    testODE.ode = problem;
    setProblem(&testODE);
    vector<double> paramList = testODE.refParam;
    // =========================

    // =========================
    // Generate a set of trajectories
    // =========================

    // REFERENCE SOLUTION
    std::vector<double> timeStepVec = {0.05};
    unsigned int nTimeSteps = 1000;
    for (size_t i = 0; i < nTimeSteps; i++) {
        timeStepVec.push_back(timeStepVec.back() * 0.999);
    }
    double finalTime = 50.0;

    // Output file
    ofstream EEResults, EIResults, MPResults, GAResults;
    string fullPath = string(DATA_PATH) + "/invMeasEE_hTest" + ".txt";
    EEResults.open(fullPath, ios::out | ofstream::trunc);

    std::vector<VectorXd> results(nTimeSteps);

    unsigned int i;
    #pragma omp parallel for num_threads(30) private(i)
    for (i = 0; i < nTimeSteps; i++) {

        double timeStep = timeStepVec[i];

        // Initialize solver
        VectorXd solution(testODE.size);
        VectorXd oldSolution(testODE.size);
        EulerBackwards EESolver(testODE, paramList);
        double time, hLoc;

        printf("timestep : %f\n", timeStep);

        int nSteps = int(finalTime / timeStep);
        oldSolution = testODE.initialCond;

        // Integrate up to final time
        time = 0.0;
        while (time < finalTime) {
            solution = EESolver.oneStep(oldSolution, timeStep);
            time = time + timeStep;
        }

        // Check if we got to the right time
        if (time > timeStep) {
            hLoc = finalTime - (time - timeStep);
            solution = EESolver.oneStep(oldSolution, hLoc);
        }

        results[i] = solution;
    }

    for (auto it : results) {
        EEResults << fixed << setprecision(50) << it.transpose() << endl;
    }

    EEResults.close();

    return 0;
}