#include <RungeKuttaSolver.hpp>
#include <iomanip>
#include <fstream>
#include <random>

int main(int argc, char* argv[])
{
    if (argc < 6) {
        throw std::invalid_argument("\nInputs:\noutput fileName \n"
                                    "Time step\n"
                                    "Final time\n"
                                    "Number of observations\n"
                                    "Noise std. dev.\n");
    }

    odeDef ODE;
    ODE.ode = FITZNAG;
    setProblem(&ODE);

    std::vector<double> param = ODE.refParam;

    Butcher tableau(RK4, EXPLICIT, 0);
    RungeKutta Method(ODE, param, tableau);

    double h = std::atof(argv[2]);
    double T = std::atof(argv[3]);
    int nObs = std::atoi(argv[4]);
    double noise = std::atof(argv[5]);

    std::default_random_engine generator{(unsigned int) time(NULL)};
    std::normal_distribution<double> noiseDist;

    int N = static_cast<int>(T / h);
    int NObsLaps = N / nObs;

    // Output file
    std::string outputFilename = std::string(DATA_PATH) + "MCMC/" + argv[1] + ".txt";
    std::ofstream output(outputFilename, std::ofstream::out | std::ofstream::trunc);

    output << nObs << std::endl;

    VectorXd solution = ODE.initialCond;
    VectorXd randomCont(solution.size());

    for (int i = 0; i < N + 1; i++) {
        solution = Method.oneStep(solution, h);
        if ((i % NObsLaps == 0) && (i > 0)) {

            for (int j = 0; j < solution.size(); j++) {
                randomCont(j) = noise * noiseDist(generator);
            }

            output << std::fixed << std::setprecision(2)
                   << h * i << "\t"
                   << std::setprecision(20)
                   << (solution + randomCont).transpose() << std::endl;
        }
    }

    output << noise << std::endl;

    output.close();

    return 0;
}
