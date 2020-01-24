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

double hamKepler(VectorXd& v, double h = 0)
{
    return 0.5 * (v(0) * v(0) + v(1) * v(1)) - 1.0 / (sqrt(v(2) * v(2) + v(3) * v(3)));
}

double hamPendulum(VectorXd& v, double h = 0)
{
    return v(0) * v(0) / 2.0 - cos(v(1));
}

double hamPendulumModif(VectorXd& v, double h)
{
    return hamPendulum(v) + h * h / 48.0 * (cos(2*v(1)) - 2*v(0)*v(0)*cos(v(1)));
}

double RPendulum(VectorXd& v, double h = 0)
{
    return (cos(2*v(1)) - 2*v(0)*v(0)*cos(v(1))) / 48.0;
}


int main(int argc, char* argv[])
{
    GetPot parser(argc, argv);

    double h = 0.1,
           T = 1000,
           p = 2.5;
    bool allTrajectory = false, RTS = true;
    int M = 10;

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
    if (parser.search("-all"))
        allTrajectory = true;
    if (parser.search("-det"))
        RTS = false;

    // ODE
    odeDef ODE;
    ODE.ode = PENDULUM;
    setProblem(&ODE);

    // Numerical method
    //Butcher tableau(STORMVER, SEPAR);
    Butcher tableau(IMPMID, IMPLICIT);
    // Butcher tableau(EXPTRAPEZ, EXPLICIT, 0);

    // Integration parameters
    auto N = static_cast<unsigned int>(T / h);

    // Random seed
    std::default_random_engine generator{(unsigned int) time(NULL)};

    std::string outputFileName = std::string(DATA_PATH) + fileName + ".txt";
    std::ofstream output(outputFileName, std::ofstream::out | std::ofstream::trunc);

    // Initialization
    RungeKuttaRandomH Method(&generator, ODE, tableau, h, p);
    RungeKutta detMethod(ODE, tableau);

    std::vector<VectorXd> solution(M, ODE.initialCond);
    VectorXd detSolution = ODE.initialCond;
    double (*hamiltonian) (VectorXd&, double) = &hamPendulum;

    if (RTS) {
        if (allTrajectory) {
            for (unsigned int i = 0; i < N; i++) {
                if (i % 5 == 0) {
                    output << std::fixed << std::setprecision(20) << h * i << "\t" << hamiltonian(detSolution, h) << "\t";
                }
                detSolution = detMethod.oneStep(h, detSolution, ODE.refParam);

                for (int j = 0; j < M; j++) {
                    if (i % 5 == 0) {
                        output << std::fixed << std::setprecision(20) << hamiltonian(solution[j], h) << "\t";
                    }
                    solution[j] = Method.oneStep(solution[j], ODE.refParam);
                }
                if (i % 5 == 0) {
                    output << std::endl;
                }
            }
        } else {
            for (unsigned int i = 0; i < N; i++) {
                detSolution = detMethod.oneStep(h, detSolution, ODE.refParam);

                for (int j = 0; j < M; j++) {
                    solution[j] = Method.oneStep(solution[j], ODE.refParam);
                }
            }
            for (int j = 0; j < M; j++) {
                output << std::fixed << std::setprecision(20) << hamiltonian(solution[j], h) << "\t";
            }
        }
    } else {
        if (allTrajectory) {
            for (unsigned int i = 0; i < N; i++) {
                if (i % 200 == 0) {
                    output << std::fixed << std::setprecision(20) << h * i << "\t" << hamiltonian(detSolution, h) << std::endl;
                }
                detSolution = detMethod.oneStep(h, detSolution, ODE.refParam);
            }
        } else {
            for (unsigned int i = 0; i < N; i++) {
                detSolution = detMethod.oneStep(h, detSolution, ODE.refParam);
            }
            output << std::fixed << std::setprecision(20) << hamiltonian(detSolution, h) << "\t";
        }
    }

    output.close();
    return 0;
}
