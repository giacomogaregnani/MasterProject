#include <RandomTimeStep.hpp>
#include <fstream>
#include <iomanip>
#include <GetPot>

double hamHenHeil(VectorXd& v, double h = 0)
{
    return 0.5 * (v(0) * v(0) + v(1) * v(1)) +
           0.5 * (v(2) * v(2) + v(3) * v(3)) +
           v(2) * v(2) * v(3) - v(3) * v(3) * v(3) / 3.0;
}

double angMomKepler(VectorXd& v, double h = 0)
{
    return v(2) * v(1) - v(3) * v(0);
}

double hamKepler(VectorXd& v, double h = 0)
{
    return 0.5 * (v(0) * v(0) + v(1) * v(1)) - 1.0 / (sqrt(v(2) * v(2) + v(3) * v(3)));
}

double hamPendulum(VectorXd& v, double h = 0)
{
    return v(0) * v(0) / 2.0 - cos(v(1));
}

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
    ODE.ode = PENDULUM;
    setProblem(&ODE);

    // Numerical method
    // Butcher tableau(GAUSS4, IMPLICIT);
    // Butcher tableau(STORMVER, SEPAR);
    Butcher tableau(IMPMID, IMPLICIT);
    // Butcher tableau(RKC, STABEXP, 150);
    // Integration parameters
    auto N = static_cast<unsigned int>(T / h);

    // Random seed
    std::default_random_engine generator{(unsigned int) time(NULL)};

    std::string outputFileName = std::string(DATA_PATH) + outputFile + ".txt";
    std::string outputDetName = std::string(DATA_PATH) + outputFile + "det.txt";
    std::ofstream outputDet(outputDetName, std::ofstream::out | std::ofstream::trunc);
    std::ofstream output(outputFileName, std::ofstream::out | std::ofstream::trunc);


    double refHam = hamPendulum(ODE.initialCond);

    // Deterministic
    RungeKutta detMethod(ODE, tableau);
    VectorXd solution = ODE.initialCond;

    // outputDet << 0.0 << "\t" << solution.transpose() << "\t" << hamPendulum(solution) << std::endl;
    outputDet << 0.0  << "\t" << 0.0 << std::endl;
    for (unsigned int i = 0; i < N; i++) {
        solution = detMethod.oneStep(h, solution, ODE.refParam);
        if (i % 10 == 0 || i * h < 10) {
            outputDet << std::fixed << std::setprecision(20)
                      // << h * (i + 1) << "\t" << solution.transpose() << "\t" << hamPendulum(solution) << std::endl;
                      << h * (i + 1) << "\t" << std::abs(hamPendulum(solution) - refHam) << std::endl;
        }
    }
    outputDet.close();

    if (chosenMethod == 0) {
        RungeKuttaRandomH Method(&generator, ODE, tableau, h, p);

        /* for (unsigned int j = 0; j < nMC; j++) {
            solution = ODE.initialCond;
            output << 0.0 << "\t" << solution.transpose() << "\t" << hamPendulum(solution) << std::endl;
            for (unsigned int i = 0; i < N; i++) {
                solution = Method.oneStep(solution, ODE.refParam);
                if (i % 100 == 0) {
                    output << std::fixed << std::setprecision(20)
                           // << h * (i + 1) << "\t" << solution.transpose() << "\t" << hamPendulum(solution) << std::endl;
                           << h * (i + 1) << "\t" << hamPendulum(solution);
                }
            }
        } */

        output << 0.0 << "\t" << 0.0 << std::endl;

        std::vector<VectorXd> solutions(nMC, ODE.initialCond);

        for (unsigned int i = 0; i < N; i++) {
            for (unsigned int j = 0; j < nMC; j++){
                solutions[j] = Method.oneStep(solutions[j], ODE.refParam);
            }

            if (i % 10 == 0 || i * h < 10) {
                double err = 0.0;
                for (unsigned int j = 0; j < nMC; j++) {
                    err += std::abs(refHam - hamPendulum(solutions[j]));
                }
                output << h * (i + 1) << "\t" << err / nMC << std::endl;
            }
        }

    } else {
        RungeKuttaAddNoise Method(&generator, ODE, tableau, h, p);

        for (unsigned int j = 0; j < nMC; j++) {
            solution = ODE.initialCond;
            output << 0.0 << "\t" << solution.transpose() << "\t" << hamPendulum(solution) << std::endl;
            for (unsigned int i = 0; i < N; i++) {
                solution = Method.oneStep(solution, ODE.refParam);
                if (i % 100 == 0) {
                    output << std::fixed << std::setprecision(20)
                           << h * (i + 1) << "\t" << solution.transpose() << "\t" << hamPendulum(solution) << std::endl;
                }
            }
        }
    }
    output.close();

    return 0;
}