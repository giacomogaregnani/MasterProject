#include <Solver.hpp>
#include "problems.hpp"
#include <fstream>
#include <iomanip>

int main(int argc, char* argv[]) {
    // =========================
    // Initialization
    // =========================
    problems problem = FITZNAG;
    odeDef testODE;
    testODE.ode = problem;
    setProblem(&testODE);
    double finalTime = 10.0;

    // =========================
    // Solution
    // =========================

    // Time step and number of trajectories
    std::vector<double> hVec = {0.1};
    for (size_t i = 0; i < 8; i++) {
        hVec.push_back(hVec.back() / 2);
    }
    std::vector<unsigned long int> NSteps;
    NSteps.push_back(100);
    for (size_t i = 0; i < 8; i++) {
        NSteps.push_back(NSteps.back() * 2);
    }
    unsigned long int M = 1000;

    // Only for implicit methods
    Butcher butcher(GAUSS, 2);
    MatrixXd A = butcher.getA();
    VectorXd b = butcher.getB();

    // Random generator
    std::default_random_engine generator{(unsigned int) time(NULL)};

    // Output file
    std::ofstream results;
    std::string finalPath = DATA_PATH + std::string("/WeakConvergence.txt");
    results.open(finalPath , std::ofstream::out | std::ofstream::trunc);

    // Reference no-noise solution
    double hRef = hVec.back() / 50;
    std::cout << "Reference solution time step : " << hRef << std::endl;
    RungeKutta SolverRef(testODE, testODE.refParam);
    VectorXd refSolution = testODE.initialCond;
    double t = 0;
    while (t < finalTime) {
        refSolution = SolverRef.oneStep(refSolution, hRef);
        t = t + hRef;
    }

    // ProbMethod realizations

    std::vector<VectorXd> MCRealizations(M, VectorXd(testODE.size));
    for (size_t k = 0; k < hVec.size(); k++) {

        double h = hVec[k];
        unsigned long int N = NSteps[k];
        std::cout << "Executing time step : " << h << std::endl;

        // Compute realizations
        unsigned long int i;
        #pragma omp parallel for num_threads(2) private(i, t)
        for (i = 0; i < M; i++) {
            ProbMethod<EulerForward> Solver(testODE, h, testODE.refParam, 0.5);
            t = 0.0;
            for (unsigned long int j = 0; j < N; j++){
                Solver.oneStep(generator, h);
                t = t + h;
            }
            printf("%f\n", t);
            MCRealizations[i] = Solver.getSolution();
        }

        // Compute Monte Carlo mean
        VectorXd MCMean = VectorXd::Zero(testODE.size);
        for (auto it : MCRealizations) {
            MCMean += it;
        }
        MCMean /= M;

        // Compute weak error and write on file
        double weakError = (refSolution - MCMean).norm();
        results << std::fixed << std::setprecision(30) << h << " " << weakError << std::endl;
    }

    results.close();

    return 0;
}
