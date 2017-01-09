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
    double hRef = 0.001;
    double finalTime = 200.0;
    int nSteps = int(finalTime / hRef);

    // random generator
    default_random_engine generator{(unsigned int) time(NULL)};

    // Write all results
    ofstream EEResults, EIResults, MPResults, GAResults;
    string fullPath = string(DATA_PATH) + "/invMeasMP_hTestSto.txt";
    EEResults.open(fullPath, ios::out | ofstream::trunc);

    // Number of realizations per solver
    size_t nRealizations = 10000;

    vector<VectorXd> fTimeEE(nRealizations);

    size_t index;
    #pragma omp parallel for num_threads(30) private(index)
    for (index = 0; index < nRealizations; index++) {
        // Initialize the solvers
        ProbMethod<MidPoint> EESolver(testODE, hRef, paramList, 0.5);

        printf("traj %zu\n", index);
        for (int i = 0; i < nSteps; i++) {
            EESolver.oneStep(generator, hRef);
        }
        fTimeEE[index] = EESolver.getSolution();
    }

    for (size_t j = 0; j < nRealizations; j++) {
        EEResults << fixed << setprecision(50) << fTimeEE[j].transpose() << endl;
    }

    EEResults.close();

    return 0;
}