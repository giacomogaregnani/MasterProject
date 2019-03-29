#include <RandomTimeStep.hpp>
#include <fstream>
#include <iomanip>
#include <GetPot>

int main(int argc, char* argv[])
{
    GetPot parser(argc, argv);

    double h = 0.1,
           T = 50,
           sigma = 0.1,
           p = 2;
    int nMC = 10;

    std::string outputFile;

    if (parser.search("-h"))
        h = parser.next(h);
    if (parser.search("-T"))
        T = parser.next(T);
    if (parser.search("-output"))
        outputFile = parser.next(" ");
    if (parser.search("-nMC"))
        nMC = parser.next(nMC);
    if (parser.search("-p"))
        p = parser.next(p);
    if (parser.search("-sigma"))
        sigma = parser.next(sigma);

    // ODE
    odeDef ODE;
    ODE.ode = LORENZ;
    setProblem(&ODE);

    // Numerical method
    Butcher tableau(EXPTRAPEZ, EXPLICIT, 0);

    // Integration parameters
    auto N = static_cast<unsigned int>(T / h);

    // Random seed
    std::default_random_engine generator{(unsigned int) time(nullptr)};

    std::string outputFileName = std::string(DATA_PATH) + outputFile + ".txt";
    std::string outputDetName = std::string(DATA_PATH) + "resultsPertICLorenz/" + outputFile + ".txt";
    std::ofstream outputDet(outputDetName, std::ofstream::out | std::ofstream::trunc);
    // std::ofstream output(outputFileName, std::ofstream::out | std::ofstream::trunc);

    // Deterministic
    RungeKutta detMethod(ODE, tableau);
    VectorXd solution = ODE.initialCond;

    double hSemiProb;
    std::uniform_real_distribution<double> unifDistribution(h - std::pow(h, p),
                                                            h + std::pow(h, p));
    std::normal_distribution<double> normalDistribution(0.0, 1.0);

    for (unsigned int j = 0; j < nMC; j++) {
        solution = ODE.initialCond;
        solution(1) += sigma * normalDistribution(generator);

        for (unsigned int i = 0; i < N; i++) {
            solution = detMethod.oneStep(h, solution, ODE.refParam);
            outputDet << h * i << "\t" << solution.transpose() << std::endl;
        }
    }
    outputDet.close();

    // Probabilistic
    /* RungeKuttaRandomH Method(&generator, ODE, tableau, h, 2.5);

    for (unsigned int j = 0; j < nMC; j++) {
        solution = ODE.initialCond;
        for (unsigned int i = 0; i < N; i++) {
            solution = Method.oneStep(solution, ODE.refParam);
            output << h * i << "\t" << solution.transpose() << std::endl;
        }
    }

    output.close(); */

    return 0;
}
