#include <RandomTimeStep.hpp>
#include <UtilitiesGG.hpp>

#include <iomanip>

double phi(VectorXd x)
{
    return x.dot(x);
}

int main(int argc, char* argv[]) {
    // Set problem and solver features
    odeDef ODE;
    ODE.ode = LORENZ;
    setProblem(&ODE);
    std::vector<double> param = ODE.refParam;
    Butcher tableau(RK4, EXPLICIT, 0);
    double T = 100, h = 1e-2;
    size_t N = static_cast<size_t>(T / h);

    std::cout << "========= TIME AVERAGE COMPUTATION (det) ============" << std::endl;

    // Open output file
    std::string fileNameTime = DATA_PATH + std::string("SolutionTimePlot.txt");
    std::ofstream solDetFile(fileNameTime, std::ofstream::out | std::ofstream::trunc);

    // Initialize structures for time integration
    double timeAverage = 0.0;
    RungeKutta detSolver(ODE, ODE.refParam, tableau);
    VectorXd detSolution = ODE.initialCond;

    // Loop over time
    for (size_t i = 0; i < N; i++) {

        detSolution = detSolver.oneStep(detSolution, h);

        timeAverage = static_cast<double>(i) / (i + 1) * timeAverage + phi(detSolution) / (i + 1);
        // solDetFile << timeAverage << "\n";
        solDetFile << detSolution.transpose() << "\n";

        if (std::fmod(i * h, T / 100.0) < h / 10.0) {
            std::cout << std::fixed << std::setprecision(1)
                      << "Time = " << i * h
                      << "s, out of " << T << "s" << std::endl;
        }
    }

    solDetFile.close();

    return 0;
}