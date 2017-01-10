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

    std::vector<double> timeStepVec = {0.1};
    unsigned int nTimeSteps = 600;
    for (size_t i = 0; i < nTimeSteps; i++) {
        timeStepVec.push_back(timeStepVec.back() * 0.999);
    }
    double finalTime = 200.0;

    // Output file
    ofstream EEResults;
    string fullPath = string(DATA_PATH) + "/invMeasEI_hTest" + ".txt";
    EEResults.open(fullPath, ios::out | ofstream::trunc);

    std::vector<VectorXd> results(nTimeSteps);

    Butcher butcherTable(RADAU, 2);
    MatrixXd A = butcherTable.getA();
    VectorXd b = butcherTable.getB();

    unsigned int i;
    #pragma omp parallel for num_threads(30) private(i)
    for (i = 0; i < nTimeSteps; i++) {

        double timeStep = timeStepVec[i];

        /*ofstream fullResults;
        string fullSolPath = string(DATA_PATH) +
                "/fullTrajLorenz_h_" +
                std::to_string(static_cast<int>(timeStep * 1e8)) +
                ".txt";
        fullResults.open(fullSolPath, ios::out | ofstream::trunc);*/

        // Initialize solver
        VectorXd solution(testODE.size);
        VectorXd oldSolution(testODE.size);
        ImplicitRK EESolver(testODE, paramList, A, b, 4);
        double time, hLoc;

        printf("timestep : %f\n", timeStep);

        oldSolution = testODE.initialCond;

        // Integrate up to final time
        time = 0.0;
        while (time < finalTime) {
            solution = EESolver.oneStep(oldSolution, timeStep);
            time = time + timeStep;
            // fullResults << oldSolution.transpose() << std::endl;
            oldSolution = solution;
        }
        // fullResults << oldSolution.transpose() << std::endl;

        // Check if we got to the right time
        if (time > timeStep) {
            hLoc = finalTime - (time - timeStep);
            solution = EESolver.oneStep(oldSolution, hLoc);
            time = time - timeStep + hLoc;
        }

        printf("nIter : %d\n", i);

        results[i] = solution;

        // fullResults.close();
    }

    for (auto it : results) {
        EEResults << fixed << setprecision(50) << it.transpose() << endl;
    }

    EEResults.close();

    return 0;
}