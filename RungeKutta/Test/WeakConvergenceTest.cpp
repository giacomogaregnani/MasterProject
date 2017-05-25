#include <RandomTimeStep.hpp>
#include <iostream>
#include <fstream>
#include <iomanip>

int main(int argc, char* argv[])
{
    odeDef ODE;
    ODE.ode = FITZNAG;
    setProblem(&ODE);

    std::vector<double> param = {0.2, 0.2, 3.0};

    Butcher tableau(EXPTRAPEZ, EXPLICIT, 0);

    double T = 1.0, h = 0.1;
    int nMC = 1e6, nExp = 6;

    std::default_random_engine generator{(unsigned int) time(NULL)};

    // Output file
    std::string outputFilename = std::string(DATA_PATH) + argv[1] + ".txt";
    std::ofstream output;
    output.open(outputFilename, std::ofstream::out | std::ofstream::trunc);

    // Reference solution
    std::string refFilename = std::string(DATA_PATH) + argv[2] + ".txt";
    std::fstream reference;
    reference.open(refFilename, std::ofstream::in);
    VectorXd refSolution(ODE.size);
    for (int i = 0; i < ODE.size; i++) {
        reference >> refSolution(i);
    }

    double exactValue = exp(param[0] * T);
    std::vector<double> errors(static_cast<unsigned int>(nExp));
    std::vector<double> orders(static_cast<unsigned int>(nExp - 1));

    // "Timestep" loop
    for (int j = 0; j < nExp; j++) {

        VectorXd MCMean = VectorXd::Zero(ODE.size);
        int N = static_cast<int>(T / h);

        std::cout << "Mean timestep = "
                  << h
                  << "\t\t"
                  << "Number of steps = "
                  << N
                  << std::endl;

        // MC loop
        int k;
        std::vector<VectorXd> MCResults(nMC);
        std::vector<double> finalTime(nMC);

        #pragma omp parallel for num_threads(30) private(k)
        for (k = 0; k < nMC; k++) {

            RungeKuttaRandomH Method(&generator, ODE, param, tableau, h, 1.0);
            VectorXd solution = ODE.initialCond;

            // Time integration loop
            for (int i = 0; i < N; i++) {
                solution = Method.oneStep(solution);
            }

            MCResults[k] = solution;
            finalTime[k] = Method.getCurrentTime();
        }

        for (auto it : MCResults) {
            MCMean += it;
        }

        double meanFinalTime = 0.0;
        for (auto it : finalTime) {
            meanFinalTime += it;
        }
        meanFinalTime /= nMC;
        std::cout << "Mean final time = " << meanFinalTime << std::endl;


        MCMean /= nMC;


        if (ODE.ode == TEST1D) {
            errors[j] = std::abs(MCMean(0) - exactValue);
            output << std::fixed << std::setprecision(20)
                   << MCMean(0) << "\t"
                   << exactValue << "\n";
        } else {
            errors[j] = std::abs(refSolution.dot(refSolution) - MCMean.dot(MCMean));
            output << std::fixed << std::setprecision(20)
                   << MCMean.dot(MCMean) << "\t"
                   << refSolution.dot(refSolution) << "\n";
        }

        h /= 2;
    }

    double meanOrder = 0.0;
    for (int j = 0; j < nExp; j++) {

        if (j > 0) {
            std::cout << "error = " << errors[j] << "\t";
            std::cout << "order = " << std::log(errors[j - 1] / errors[j]) / std::log(2.0) << std::endl;
            meanOrder += std::log(errors[j - 1] / errors[j]) / std::log(2.0);
        } else {
            std::cout << "error = " << errors[j] << std::endl;
        }
    }
    meanOrder /= nExp - 1;
    std::cout << "Mean order of convergence = " << meanOrder << std::endl;

    output.close();

    return 0;
}
