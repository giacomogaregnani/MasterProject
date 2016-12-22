#include <Solver.hpp>
#include "problems.hpp"
#include <fstream>
#include <cmath>

int main(int argc, char* argv[]) {
    // =========================
    // INITIALIZATION
    // =========================
    problems problem = LORENZ;

    // Set problem
    odeDef testODE;
    testODE.ode = problem;
    setProblem(&testODE);
    //

    std::vector<double> paramList = testODE.refParam;
    double finalTime = 100.0;
    // =========================

    // =========================
    // SOLVE USING GAUSS
    // =========================

    double h = 0.01;
    int M = 0;

    Butcher butcher(GAUSS, 2);
    MatrixXd A = butcher.getA();
    VectorXd b = butcher.getB();

    std::default_random_engine generator{(unsigned int) time(NULL)};

    std::ofstream results;
    std::string finalPath = DATA_PATH + std::string("/LorenzTrajSetLongTime.txt");
    results.open(finalPath , std::ofstream::out | std::ofstream::trunc);

    // Reference no-noise solution
    ImplicitRK SolverRef(testODE, testODE.refParam, A, b, 4);
    VectorXd solution = testODE.initialCond;
    double t = 0;
    while (t < finalTime) {
        solution = SolverRef.oneStep(solution, h);
        t = t + h;
        results << t << " " << solution.transpose() << std::endl;
    }

    // ProbMethod realizations
    impProbMethod Solver(testODE, h, testODE.refParam, 0.5, A, b, 4);
    for (int i = 0; i < M; i++) {
        t = 0;
        while (t < finalTime) {
            Solver.oneStep(generator, h);
            t = t + h;
            results << t << " " << Solver.getSolution().transpose() << std::endl;
        }
        Solver.resetIC();
    }

    results.close();

    return 0;
}
