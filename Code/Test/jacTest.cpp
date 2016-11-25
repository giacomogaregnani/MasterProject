#include <fstream>
#include "problems.hpp"

int main(int argc, char* argv[]) {
    // =========================
    // INITIALIZATION
    // =========================
    problems problem = HIRES;
    // Set problem
    odeDef testODE;
    testODE.ode = problem;
    setProblem(&testODE);
    // =========================

    // Write all results
    std::string folder(DATA_PATH);
    std::string slash("/");
    std::ofstream fullResults;
    std::string fullResPath = "jacobianHires.txt";
    std::string fullPath = folder + slash + fullResPath;
    fullResults.open(fullPath, std::ofstream::out | std::ofstream::trunc);

    MatrixXd initialVar = 1e-3 * MatrixXd::Identity(testODE.size, testODE.size);

    fullResults << testODE.odeJac(VectorXd::Random(testODE.size), testODE.refParam) << "\n";
    fullResults << initialVar;

    fullResults.close();
    return 0;
}