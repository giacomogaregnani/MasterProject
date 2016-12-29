#include <Solver.hpp>
#include "problems.hpp"
#include "mcmcTools.hpp"
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
    double finalTime = 20.0;

    // =========================
    // Solution
    // =========================

    unsigned long int M = 5;

    // Random generator
    std::default_random_engine generator{(unsigned int) time(NULL)};
    std::normal_distribution<double> normal(0.0, 1.0);

    // Output file
    std::ofstream results;
    std::string finalPath = DATA_PATH + std::string("/Adaptivity.txt");
    results.open(finalPath , std::ofstream::out | std::ofstream::trunc);

    // ProbMethod realizations

    double t = 0, h = 0.1, targetVar = 0.01, varNorm, sigma = 0.5;
    int odeSize = testODE.size;
    // Compute realizations
    //std::vector<ProbMethod<EulerForward>> Solver(M, ProbMethod<EulerForward>(testODE, h, testODE.refParam, 0.5));
    std::vector<VectorXd> currentSolution(M, VectorXd(testODE.size)), attemptSolution(M, VectorXd(testODE.size));
    VectorXd mean, variances, disturb(testODE.size);
    while (t < finalTime) {

        if (t + h > finalTime) {
            h = finalTime - t;
        }

        for (unsigned int i = 0; i < M; i++) {
            for (int j = 0; j < odeSize; j++) {
                disturb(j) = normal(generator);
            }
            attemptSolution[i] = currentSolution[i] +
                    h * testODE.odeFunc(currentSolution[i], testODE.refParam) +
                    pow(h, 1.5) * sigma * disturb;
        }

        mean = computeMeans(attemptSolution);
        variances = computeVariances(attemptSolution, mean);
        varNorm = variances.norm();

        if (varNorm < targetVar) {
            t = t + h;
            for (unsigned int i = 0; i < M; i++) {
                currentSolution[i] = attemptSolution[i];
            }
        }
        h = h * std::min(2.0, std::max(pow(targetVar / varNorm, 2.0), 0.5));


        // Update the time step

        /*for (auto it : currentSolution) {
            std::cout << it.transpose() << std::endl;
        }
        std::cout << " =========== " << std::endl; */
        std::cout << "time " << t << " h = " << h << " var = " << varNorm << std::endl;
        std::cout << " =========== " << std::endl;

        results << std::fixed << std::setprecision(30) << t << " " << varNorm << std::endl;
    }


    results.close();

    return 0;
}
