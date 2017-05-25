#include <RungeKuttaSolver.hpp>
#include <iostream>
#include <iomanip>
#include <fstream>

int main(int argc, char* argv[])
{
    if (argc < 4) {
        throw std::invalid_argument("Input output fileName \n"
                                    "Time step\n"
                                    "Final time\n");
    }


    odeDef ODE;
    ODE.ode = FITZNAG;
    setProblem(&ODE);

    std::vector<double> param = ODE.refParam;

    Butcher tableau(RK4, EXPLICIT, 50);
    RungeKutta Method(ODE, param, tableau);

    double T = std::atof(argv[3]);
    double h = std::atof(argv[2]);
    int N = static_cast<int>(T / h);

    // Output file
    std::string outputFilename = std::string(DATA_PATH) + argv[1] + ".txt";
    std::ofstream output(outputFilename, std::ofstream::out | std::ofstream::trunc);

    VectorXd solution = ODE.initialCond;
    for (int i = 0; i < N; i++) {
        solution = Method.oneStep(solution, h);
    }

    output << std::fixed << std::setprecision(20)
           << solution.transpose() << std::endl;
    output.close();

    return 0;
}

