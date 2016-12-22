#include <cmath>
#include <vector>
#include <random>
#include <Solver.hpp>
#include <fstream>
#include "problems.hpp"
#include <iomanip>

int main(void)
{
    // =========================
    // INITIALIZATION
    // =========================
    problems problem = LORENZ;
    odeDef testODE;
    testODE.ode = problem;
    setProblem(&testODE);
    std::vector<double> paramList = testODE.refParam;
    std::string filepath("varMCResults_BE_Lorenz");

    // =========================
    // IMPORT EXACT SOLUTION
    // =========================
    std::string refFile("refSolLorenz");
    std::fstream refSolution;
    refSolution.open(std::string(DATA_PATH) + "/" + refFile + ".txt", std::ios_base::in);
    double finalTime, unused;
    refSolution >> finalTime;
    refSolution >> unused;
    refSolution >> unused;
    VectorXd reference(testODE.size);
    for (int j = 0; j < testODE.size; j++) {
        refSolution >> reference(j);
    }
    double refVal = reference.dot(reference);

    // =========================
    // SOLUTION GENERATION
    // =========================

    // h and M
    std::vector<double> h;
    double hMax = 0.2;
    for (int i = 0; i < 10; i++) {
        h.push_back(hMax);
        hMax /= 2.0;
    }
    std::vector<unsigned int> M = {};
    unsigned int Mmin = 1;
    for (int i = 0; i < 10; i++) {
        M.push_back(Mmin);
        Mmin *= 2;
    }
    std::default_random_engine generator{(unsigned int) time(NULL)};

    // Number of experiences
    unsigned int nExp = 100;

    // temporary vector and sum
    VectorXd tmp(testODE.size);
    double sum;

    // MC estimator
    std::vector<double> MC(nExp, 0.0);

    // Output file
    std::ofstream results;
    std::string finalpath = std::string(DATA_PATH) + "/" + filepath + ".txt";
    results.open(finalpath, std::ofstream::out | std::ofstream::trunc);

    // Butcher tableau
    Butcher butcher(RADAU, 2);
    MatrixXd A = butcher.getA();
    VectorXd b = butcher.getB();

    for (auto it : h) {
        std::cout << "Executing timestep : " << it << std::endl;

        // Compute the MC estimator
        double mean = 0, tmpSq, variance = 0, meanErr = 0;

        unsigned int j;
        #pragma omp parallel for num_threads(30) private(j)
        for (j = 0; j < nExp; j++) {
            impProbMethod refSolver(testODE, it, paramList, 0.1, A, b, 4);
            MC[j] = 0.0;
            for (unsigned int i = 0; i < M.front(); i++) {
                refSolver.resetIC();
                double time = 0;
                while (time < finalTime) {
                    refSolver.oneStep(generator, it);
                    time = time + it;
                }
                tmp = refSolver.getSolution();
                MC[j] += tmp.dot(tmp);
            }
            MC[j] /= M.front();
        }

        for (auto itMean : MC) {
            std::cout << itMean << " " << refVal << std::endl;
            mean += itMean;
        }
        mean /= nExp;

        // COMPUTE BIAS SQUARED
        meanErr = mean - refVal;
        meanErr *= meanErr;

        // Compute variance of the MC estimator
        for (j = 0; j < nExp; j++) {
            tmpSq = MC[j] - mean;
            variance += tmpSq * tmpSq;
        }
        variance /= nExp;

        // Print result on file
        results << it << " " << std::fixed << std::setprecision(30) << variance << " " << meanErr << "\n";
    }

    results.close();

    return 0;
}


