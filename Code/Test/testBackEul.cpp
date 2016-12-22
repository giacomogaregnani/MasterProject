#include <Solver.hpp>
#include "problems.hpp"
#include <fstream>
#include <cmath>

int main(int argc, char* argv[]) {
    // =========================
    // INITIALIZATION
    // =========================
    problems problem = POISSON;

    // Set problem
    odeDef testODE;
    testODE.ode = problem;
    setProblem(&testODE);
    //

    std::vector<double> paramList = testODE.refParam;
    double finalTime = 1.0;
    // =========================

    // =========================
    // REFERENCE USING RK4
    // =========================

    double h = 0.000001;
    double time = 0;

    VectorXd refSol = testODE.initialCond;

    RungeKutta refSolver(testODE, testODE.refParam);

    while (time < finalTime) {
        refSol = refSolver.oneStep(refSol, h);
        time = time + h;
    }

    // =========================
    // SOLVE USING BACKWARD EULER
    // =========================

    double hBE = 0.5;
    std::vector<double> error;

    Butcher butcher(RADAU, 2);
    MatrixXd A = butcher.getA();
    VectorXd b = butcher.getB();

    std::cout << A
              << std::endl
              << "========"
              << std::endl
              << b.transpose()
              << std::endl;

    ImplicitRK Solver(testODE, testODE.refParam, A, b, 4);
    for (int i = 0; i < 14; i++) {
        VectorXd solution = testODE.initialCond;
        std::cout << i << std::endl;
        time = 0;
        while (time < finalTime) {
            solution = Solver.oneStep(solution, hBE);
            time = time + hBE;
        }
        hBE /= 2;
        error.push_back((refSol - solution).norm());
    }

    std::ofstream results;
    std::string finalpath = DATA_PATH + std::string("/errorsimplicit.txt");
    results.open(finalpath , std::ofstream::out | std::ofstream::trunc);

    for (auto it : error) {
        results << it << "\n";
    }

    results.close();

    return 0;
}
