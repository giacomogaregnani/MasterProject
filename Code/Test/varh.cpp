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
    // testODE.refParam[0] = 0.5;
    vector<double> paramList = testODE.refParam;
    // =========================

    // =========================
    // Generate a set of trajectories
    // =========================

    size_t nTimeSteps = 1;
    std::vector<double> h = {0.01};
    for (size_t i = 0; i < nTimeSteps - 1; i++) {
        h.push_back(h.back() / 2.0);
    }

    double finalTime = 40.0;
    unsigned int nTrajectories = 100;

    std::default_random_engine generator{(unsigned int) time(NULL)};

    Butcher butcher(RADAU, 2);
    MatrixXd A = butcher.getA();
    VectorXd b = butcher.getB();

    for (auto meanTimeStep : h) {

        std::cout << meanTimeStep << std::endl;

        // Output file
        ofstream EEResults;
        string fullPath = string(DATA_PATH) + "/invMeasMPLorenz" +
                std::to_string(static_cast<int>(meanTimeStep * 1e5)) + ".txt";
        EEResults.open(fullPath, ios::out | ofstream::trunc);
        ofstream fullResults;
        string fullFullPath = string(DATA_PATH) + "/invMeasMPLorenzFull" +
                          std::to_string(static_cast<int>(meanTimeStep * 1e5)) + ".txt";
        fullResults.open(fullFullPath, ios::out | ofstream::trunc);

        // Generator of time steps
        // std::uniform_real_distribution<double> hDistribution(0.0, 2 * meanTimeStep);
        std::uniform_real_distribution<double> hDistribution(meanTimeStep - pow(meanTimeStep, 2.5),
                                                             meanTimeStep + pow(meanTimeStep, 2.5));


        // All trajectories
        std::vector<VectorXd> results(nTrajectories);

        int meanNSteps = static_cast<int> (finalTime / meanTimeStep);

        unsigned int j = 0;
        //#pragma omp parallel for num_threads(30) private(j)
        for (j = 0; j < nTrajectories; j++) {

            // Initialize solver
            VectorXd solution(testODE.size);
            VectorXd oldSolution(testODE.size);
            VectorXd oldOldSolution(testODE.size);
            double time, hLoc;
            MidPoint EESolver(testODE, paramList);
            oldSolution = testODE.initialCond;

            /* hLoc = hDistribution(generator);
            solution = EESolver.oneStep(oldSolution, hLoc); */

            // Integrate up to final time
             time = 0.0;

            fullResults << "0" << "\t" << testODE.initialCond.transpose() << std::endl;

            for (int i = 0; i < meanNSteps; i++) {
            // while (time < finalTime) {
                hLoc = hDistribution(generator);
                solution = EESolver.oneStep(oldSolution, hLoc);
                time = time + hLoc;
                oldOldSolution = oldSolution;
                fullResults << time << "\t" << solution.transpose() << std::endl;
                oldSolution = solution;
            }

            // Check if we got to the right time
            /* if (time > finalTime) {
                hLoc = finalTime - (time - hLoc);
                solution = EESolver.oneStep(oldOldSolution, hLoc);
            } */

            results[j] = solution;
        }

        for (auto iterator : results) {
            EEResults << std::fixed << std::setprecision(30) << iterator.transpose() << std::endl;
        }

        EEResults.close();
        fullResults.close();
    }

    return 0;
}