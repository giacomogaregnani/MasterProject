#include <RandomTimeStep.hpp>
#include <iostream>
#include <fstream>
#include <iomanip>

int main(int argc, char* argv[])
{
    odeDef ODE;
    ODE.ode = TEST1D;
    setProblem(&ODE);

    std::vector<double> param = {-1.0};

    Butcher tableau(EXPTRAPEZ, EXPLICIT, 0);

    double T = 1;
    double h = 0.1;
    int N, nExp = 8, nMC = 10;
    std::vector<double> error(nExp);

    std::default_random_engine generator{(unsigned int) time(NULL)};

    double exactValue = exp(param[0] * T);

    // Output file
    std::string outputFilename = std::string(DATA_PATH) + argv[1] + ".txt";
    std::ofstream output;

    output.open(outputFilename, std::ofstream::out | std::ofstream::trunc);

    // "Timestep" loop
    for (int j = 0; j < nExp; j++) {

        std::cout << "Mean timestep = " << h << std::endl;

        RungeKuttaRandomH Method(&generator, ODE, param, tableau, h, 2.5);
        double MCError = 0.0, finalTime = 0.0;

        // MC loop
        for (int k = 0; k < nMC; k++) {
            VectorXd solution = ODE.initialCond;
            N = static_cast<int>(T / h);

            // Time integration loop
            for (int i = 0; i < N; i++) {
                solution = Method.oneStep(solution);
            }

            finalTime += Method.getCurrentTime();
            MCError += std::abs(solution(0) - exactValue);
            Method.resetTimeToZero();
        }

        MCError /= nMC;
        finalTime /= nMC;
        std::cout << "Mean final time = " << finalTime << std::endl;
        error[j] = MCError;
        h /= 2;

        output << std::fixed << std::setprecision(20) << h << "\t" << MCError << "\n";


    }

    for (int j = 0; j < nExp - 1; j++) {
        std::cout << "order = " << std::log(error[j] / error[j + 1]) / std::log(2.0) << std::endl;
    }

    return 0;
}