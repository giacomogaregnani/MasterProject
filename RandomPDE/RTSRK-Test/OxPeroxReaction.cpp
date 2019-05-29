#include <RandomTimeStep.hpp>
#include <fstream>
#include <iomanip>
#include <GetPot>

int main(int argc, char* argv[])
{
    GetPot parser(argc, argv);

    double h = 0.1, p = 1.5,
            T = 50;
    int nMC = 10;
    int chosenMethod = 0;

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
    if (parser.search("-method"))
        chosenMethod = parser.next(chosenMethod);

    // ODE
    odeDef ODE;
    ODE.ode = PEROX;
    setProblem(&ODE);

    // Numerical method
    Butcher tableau(RKC, STABEXP, 150);
    // Integration parameters
    auto N = static_cast<unsigned int>(T / h);

    // Random seed
    std::default_random_engine generator{(unsigned int) time(NULL)};

    std::string outputFileName = std::string(DATA_PATH) + outputFile + ".txt";
    std::ofstream output(outputFileName, std::ofstream::out | std::ofstream::trunc);

    // Deterministic
    RungeKutta detMethod(ODE, tableau);
    VectorXd solution = ODE.initialCond;

    if (chosenMethod == 0) {
        RungeKuttaRandomH Method(&generator, ODE, tableau, h, p);
        for (unsigned int j = 0; j < nMC; j++) {
            solution = ODE.initialCond;
            output << 0.0 << "\t" << solution.transpose() << std::endl;
            for (unsigned int i = 0; i < N; i++) {
                solution = Method.oneStep(solution, ODE.refParam);
                if (i % 10 == 0) {
                    output << std::fixed << std::setprecision(20)
                           << h * (i + 1) << "\t" << solution.transpose() <<  std::endl;
                }
            }
        }
    } else {
        RungeKuttaAddNoise Method(&generator, ODE, tableau, h, p, 0.0001);
        for (unsigned int j = 0; j < nMC; j++) {
            solution = ODE.initialCond;
            output << 0.0 << "\t" << solution.transpose() << std::endl;
            for (unsigned int i = 0; i < N; i++) {
                solution = Method.oneStep(solution, ODE.refParam);
                if (i % 10 == 0) {
                    output << std::fixed << std::setprecision(20)
                           << h * (i + 1) << "\t" << solution.transpose() <<  std::endl;
                }
            }
        }
    }
    output.close();

    return 0;
}