#include <RungeKuttaSolver.hpp>
#include <iomanip>
#include <fstream>
#include <random>
#include <GetPot>

int main(int argc, char* argv[])
{
    // Parse input parameters
    GetPot parser(argc, argv);
    std::string obsFileName = "observations";
    double h = 0.001;
    double T = 1;
    int N = 100;
    double theta = 0.01, noise = 0.01;

    if (parser.search("-obsFile"))
        obsFileName = parser.next("observations");
    if (parser.search("-h"))
        h = parser.next(h);
    if (parser.search("-T"))
        T = parser.next(T);
    if (parser.search("-N"))
        N = parser.next(N);
    if (parser.search("-theta"))
        theta = parser.next(theta);
    if (parser.search("-noise"))
        noise = parser.next(noise);

    odeDef ODE;
    ODE.ode = HEAT;
    setProblem(&ODE, N);

    std::vector<double> param = {theta};

    Butcher tableau(IMPEULER, IMPLICIT, 0);
    RungeKutta Method(ODE, param, tableau);

    std::default_random_engine generator{(unsigned int) time(NULL)};
    std::normal_distribution<double> noiseDist;

    int nTimeSteps = static_cast<int>(T / h);

    // Output file
    std::string outputFilename = std::string(DATA_PATH) + "HeatEquation/" + obsFileName + ".txt";
    std::ofstream output(outputFilename, std::ofstream::out | std::ofstream::trunc);

    // Build initial condition (object of inference)
    double dX = 1.0 / N;

    VectorXd initialCond = VectorXd::Zero(ODE.size);
    for (int i = 0; i < ODE.size; i++) {
        if ((i + 1) * dX < 0.9 && (i + 1) * dX > 0.7)
            initialCond(i) = 1;
    }

    VectorXd solution = initialCond;
    VectorXd randomCont(solution.size());


    for (int i = 0; i < nTimeSteps + 1; i++) {
        solution = Method.oneStep(solution, h);
    }

    for (int i = 0; i < ODE.size; i++) {
        randomCont(i) = noise * noiseDist(generator);
    }
    solution += randomCont;

    output << noise << std::endl;
    output << N << std::endl;
    output << solution.transpose() << std::endl;

    output.close();

    return 0;
}
