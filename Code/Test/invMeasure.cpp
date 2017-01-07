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
    double hRef = 0.1;
    double finalTime = 50.0;
    int nSteps = int(finalTime / hRef);

    Butcher butcher(GAUSS, 2);
    MatrixXd A = butcher.getA();
    MatrixXd b = butcher.getB();

    // random generator
    default_random_engine generator{(unsigned int) time(NULL)};

    for (int k = 0; k < 5; k++) {

        hRef /= 2;
        cout << "h = " << hRef << endl;

        // Write all results
        ofstream EEResults, EIResults, MPResults, GAResults;
        string fullPath = string(DATA_PATH) + "/invMeasEE_h" + to_string((int) (hRef*1e6)) + ".txt";
        EEResults.open(fullPath, ios::out | ofstream::trunc);
        fullPath = string(DATA_PATH) + "/invMeasEI_h" + to_string((int) (hRef*1e6)) + ".txt";
        EIResults.open(fullPath, ios::out | ofstream::trunc);
        fullPath = string(DATA_PATH) + "/invMeasMP_h" + to_string((int) (hRef*1e6)) + ".txt";
        MPResults.open(fullPath, ios::out | ofstream::trunc);
        fullPath = string(DATA_PATH) + "/invMeasGA_h" + to_string((int) (hRef*1e6)) + ".txt";
        GAResults.open(fullPath, ios::out | ofstream::trunc);

        // Number of realizations per solver
        size_t nRealizations = 10000;

        vector<VectorXd> fTimeEE(nRealizations);
        vector<VectorXd> fTimeEI(nRealizations);
        vector<VectorXd> fTimeMP(nRealizations);
        vector<VectorXd> fTimeGA(nRealizations);

        size_t index;
        #pragma omp parallel for num_threads(30) private(index)
        for (index = 0; index < nRealizations; index++) {
            // Initialize the solvers
            ProbMethod<EulerForward> EESolver(testODE, hRef, paramList, 0.5);
            ProbMethod<EulerBackwards> EISolver(testODE, hRef, paramList, 0.5);
            ProbMethod<MidPoint> MPSolver(testODE, hRef, paramList, 0.5);
            impProbMethod GASolver(testODE, hRef, paramList, 0.5, A, b, 4);


            printf("traj %zu\n", index);
            for (int i = 0; i < nSteps; i++) {
                EESolver.oneStep(generator, hRef);
                EISolver.oneStep(generator, hRef);
                MPSolver.oneStep(generator, hRef);
                GASolver.oneStep(generator, hRef);
            }
            fTimeEE[index] = EESolver.getSolution();
            fTimeEI[index] = EISolver.getSolution();
            fTimeMP[index] = MPSolver.getSolution();
            fTimeGA[index] = GASolver.getSolution();
        }

        for (size_t j = 0; j < nRealizations; j++) {
            EEResults << fixed << setprecision(50) << fTimeEE[j].transpose() << endl;
            EIResults << fixed << setprecision(50) << fTimeEI[j].transpose() << endl;
            MPResults << fixed << setprecision(50) << fTimeMP[j].transpose() << endl;
            GAResults << fixed << setprecision(50) << fTimeGA[j].transpose() << endl;
        }

        EEResults.close();
        EIResults.close();
        MPResults.close();
        GAResults.close();

    }



    return 0;
}