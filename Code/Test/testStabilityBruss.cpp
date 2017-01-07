#include <Solver.hpp>
#include "problems.hpp"
#include <fstream>
#include <iomanip>

int main(int argc, char* argv[]) {
    // =========================
    // Initialization
    // =========================
    problems problem = BRUSS;
    odeDef testODE;
    testODE.ode = problem;
    setProblem(&testODE);
    double finalTime = 50.0;

    // =========================
    // Solution
    // =========================

    // Stiffness index
    double lambda = 4.0 * testODE.refParam[0] * ((testODE.size / 2.0 + 1) * (testODE.size / 2.0 + 1));

    // Random generator
    std::default_random_engine generator{(unsigned int) time(NULL)};

    // Output file
    std::ofstream results;
    std::string finalPath = DATA_PATH + std::string("/brussSol.txt");
    results.open(finalPath , std::ofstream::out | std::ofstream::trunc);

    double h = 0.9 / (2 * lambda);
    std::cout << "Required time step : " << h << std::endl;

    EulerForward Solver(testODE, testODE.refParam);
    unsigned long int N = static_cast<unsigned long int> (finalTime / h);
    double t = 0.0;
    VectorXd solution = testODE.initialCond;
    /*results << t << "\t" << solution.transpose() << std::endl;
    for (unsigned long int j = 0; j < N; j++){
        solution = Solver.oneStep(solution, h);
        results << t << "\t" << solution.transpose() << std::endl;
        t = t + h;
    }*/

    h = 0.001;
    int stages = static_cast<int>(std::ceil(sqrt(0.5 * h * lambda))) + 1;
    unsigned long int nStab = static_cast<unsigned long int> (finalTime / h);
    std::cout << "Required number of stages : " << stages << std::endl;
    RKC stabSolver(testODE, testODE.refParam, stages, 0.0);
    t = 0.0;
    solution = testODE.initialCond;
    results << t << "\t" << solution.transpose() << std::endl;
    for (unsigned long int j = 0; j < nStab; j++){
        solution = stabSolver.oneStep(solution, h, stages);
        results << t << "\t" << solution.transpose() << std::endl;
        t = t + h;
    }

    results.close();

    return 0;
}
