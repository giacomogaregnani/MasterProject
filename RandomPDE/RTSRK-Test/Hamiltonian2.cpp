#include <RandomTimeStep.hpp>
#include <fstream>
#include <iomanip>
#include <GetPot>

double HenHeilHamiltonian(VectorXd& v)
{
    return 0.5 * (v(0) * v(0) + v(1) * v(1)) +
           0.5 * (v(2) * v(2) + v(3) * v(3)) +
           v(2) * v(2) * v(3) - v(3) * v(3) * v(3) / 3.0;
}

double hamKepler(VectorXd& v)
{
    return 0.5 * (v(0) * v(0) + v(1) * v(1)) - 1.0 / (sqrt(v(2) * v(2) + v(3) * v(3)));
}

double hamPendulum(VectorXd& v)
{
    return v(0) * v(0) / 2.0 - cos(v(1));
}

int main(int argc, char* argv[])
{
    GetPot parser(argc, argv);

    double h = 0.1,
            T = 1000,
            p = 2.5;
    bool RTS = true;
    int M = 10, nT = 4;

    std::string fileName;

    if (parser.search("-h"))
        h = parser.next(h);
    if (parser.search("-T"))
        T = parser.next(T);
    if (parser.search("-output"))
        fileName = parser.next(" ");
    if (parser.search("-M"))
        M = parser.next(M);
    if (parser.search("-p"))
        p = parser.next(p);
    if (parser.search("-det"))
        RTS = false;
    if (parser.search("-nT"))
        nT = parser.next(nT);

    // ODE
    odeDef ODE;
    ODE.ode = PENDULUM;
    setProblem(&ODE);

    // Numerical method
    Butcher tableau(IMPMID, IMPLICIT, 0);

    // Random seed
    std::default_random_engine generator{(unsigned int) time(NULL)};

    std::string outputFileName = std::string(DATA_PATH) + fileName + ".txt";
    std::ofstream output(outputFileName, std::ofstream::out | std::ofstream::trunc);

    // Initialization
    RungeKuttaRandomH Method(&generator, ODE, tableau, h, p);
    RungeKutta detMethod(ODE, tableau);

    std::vector<VectorXd> solution(M, ODE.initialCond);
    VectorXd detSolution = ODE.initialCond;
    double (*hamiltonian) (VectorXd&) = &hamPendulum;

    T /= 2;

    if (RTS) {
        for (int i = 0; i < nT; i++) {
            T = 2 * T;
            auto N = static_cast<int>(T / h);
            std::cout << "next final time = " << T << std::endl;

            for (int j = 0; j < M; j++) {
                for (int k = 0; k < N; k++) {
                    solution[j] = Method.oneStep(solution[j], ODE.refParam);
                }
            }
            std::cout << "current time = " << T << std::endl;

            output << T << "\t";
            for (int j = 0; j < M; j++) {
                output << std::fixed << std::setprecision(20) << hamiltonian(solution[j]) << "\t";
            }
            output << std::endl;
        }
    } else {
        for (int i = 0; i < nT; i++) {
            T = 2 * T;
            auto N = static_cast<int>(T / h);
            std::cout << "next final time = " << T << std::endl;

            for (int j = 0; j < N; j++) {
                detSolution = Method.oneStep(detSolution, ODE.refParam);
            }
            std::cout << "current time = " << T << std::endl;

            output << std::fixed << std::setprecision(20) << T << "\t" << hamiltonian(detSolution) << std::endl;
        }
    }

    output.close();
    return 0;
}
