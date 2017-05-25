#include <RungeKuttaSolver.hpp>
#include <iostream>
#include <iomanip>
#include <fstream>

int main(int argc, char* argv[])
{
    if (argc < 5) {
        throw std::invalid_argument("Input output fileName \n"
                                    "Time step\n"
                                    "Final time\n"
                                    "Number of observations\n");
    }

    odeDef ODE;
    ODE.ode = LORENZ;
    setProblem(&ODE);

    std::vector<double> param = ODE.refParam;

    Butcher tableau(RK4, EXPLICIT, 0);
    RungeKutta Method(ODE, param, tableau);

    double h = std::atof(argv[2]);
    double T = std::atof(argv[3]);
    int nObs = std::atoi(argv[4]);
    int N = static_cast<int>(T / h);
    int NObsLaps = N / nObs;

    // Output file
    std::string outputFilename = std::string(DATA_PATH) + argv[1] + ".txt";
    std::ofstream output(outputFilename, std::ofstream::out | std::ofstream::trunc);

    output << nObs << std::endl;

    VectorXd solution = ODE.initialCond;
    for (int i = 0; i < N + 1; i++) {
        solution = Method.oneStep(solution, h);
        if ((i % NObsLaps == 0) && (i > 0)) {
            output << std::fixed << std::setprecision(2)
                   << h * i << "\t"
                   << std::setprecision(20)
                   << solution.transpose() << std::endl;
        }
    }

    output.close();

    return 0;
}
