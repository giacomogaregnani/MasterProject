#include <RandomTimeStep.hpp>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>

int main(int argc, char* argv[])
{
    odeDef ODE;
    ODE.ode = FITZNAG;
    setProblem(&ODE);

    std::vector<double> param = {0.2, 0.2, 3.0};

    Butcher tableau(EULERFORWARD, EXPLICIT, 0);

    double T = 10, h = 0.1;
    int N = static_cast<int>(T / h);
    int nMC = 100;

    std::default_random_engine generator{(unsigned int) time(NULL)};

    // Output file
    std::string outputFilename = std::string(DATA_PATH) + argv[1] + ".txt";
    std::ofstream output;

    output.open(outputFilename, std::ofstream::out | std::ofstream::trunc);

    // "Timestep" loop
    std::cout << "Mean timestep = " << h << std::endl;

    RungeKuttaRandomH Method(&generator, ODE, param, tableau, h, 1.5);

    // MC loop
    for (int k = 0; k < nMC; k++) {

        VectorXd solution = ODE.initialCond;

        // Time integration loop
        for (int i = 0; i < N; i++) {
            solution = Method.oneStep(solution);
            output << std::fixed << std::setprecision(20)
                   << Method.getCurrentTime()
                   << "\t"
                   << solution.transpose()
                   << "\n";
        }

        Method.resetTimeToZero();

    }

    output.close();

    return 0;
}
