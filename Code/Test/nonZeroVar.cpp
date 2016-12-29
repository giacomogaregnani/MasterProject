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
    problems problem = FITZNAG;
    odeDef testODE;
    testODE.ode = problem;
    setProblem(&testODE);
    std::vector<double> paramList = testODE.refParam;
    std::string filepath("varMCResults_nonZVar");

    // =========================
    // IMPORT EXACT SOLUTION
    // =========================
    std::string refFile("refSolFitznag");
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

    // Number of experiences
    unsigned int nExp = 10;

    // h and M
    std::vector<double> hVec = {0.2};
    for (unsigned int i = 0; i < 8; i++) {
        hVec.push_back(hVec.back() / 2.0);
    }

    std::vector<int> MVec = {10};
    for (int i = 0; i < 3; i++) {
        MVec.push_back(MVec.back() * 2);
    }

    // random seed
    std::default_random_engine generator{(unsigned int) time(NULL)};

    // temporary vector and sum
    VectorXd tmp(testODE.size);

    // MC estimator
    std::vector<double> MC(nExp, 0.0);

    // Output file
    std::ofstream results;
    std::string finalpath = std::string(DATA_PATH) + "/" + filepath + ".txt";
    results.open(finalpath, std::ofstream::out | std::ofstream::trunc);

    // save initial condition
    VectorXd initCond = testODE.initialCond;
    VectorXd rnd(testODE.size);

    for (auto h : hVec) {
        for (auto M : MVec) {
            std::cout << "Executing timestep : " << h << std::endl;

            // Compute the MC estimator
            double mean = 0, tmpSq, variance = 0, meanErr = 0;

            unsigned int j;
            for (j = 0; j < nExp; j++) {
                MC[j] = 0.0;
                for (unsigned int i = 0; i < M; i++) {
                    rnd = 0.1 * VectorXd::Random(testODE.size);
                    testODE.initialCond = initCond + rnd;
                    // std::cout << testODE.initialCond.transpose() << std::endl;
                    ProbMethod<EulerForward> refSolver(testODE, h, paramList, 0.5);
                    double time = 0;
                    while (time < finalTime) {
                        refSolver.oneStep(generator, h);
                        time = time + h;
                    }
                    tmp = refSolver.getSolution();
                    MC[j] += tmp.dot(tmp);
                }
                MC[j] /= M;
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
            results << h << " " << M << " " << std::fixed << std::setprecision(30) << variance << " " << meanErr << "\n";
        }
    }

    results.close();

    return 0;
}


