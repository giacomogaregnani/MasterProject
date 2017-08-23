#include <RandomTimeStep.hpp>
#include <iostream>
#include <fstream>
#include <iomanip>

double todaHamiltonian(VectorXd& v)
{
    return 0.5 * (v(0) * v(0) + v(1) * v(1) + v(2) * v(2))
           + exp(v(3) - v(4)) + exp(v(4) - v(5)) + exp(v(5) - v(3));
}

int main(int argc, char* argv[])
{
    if (argc < 5) {
        throw std::invalid_argument("Number of inputs must be 4: \n"
                                            "Time step\n"
                                            "Final time\n"
                                            "File name for output\n"
                                            "exponent p\n"
                                            "==========================");
    }

    // ODE
    odeDef ODE;
    ODE.ode = TODA;
    setProblem(&ODE);
    std::vector<double> param = ODE.refParam;

    // Numerical method
    Butcher tableau(GAUSS4, IMPLICIT, 0);
    Butcher tableauRef(GAUSS4, IMPLICIT, 0);
    double p = std::atof(argv[4]);

    // Integration parameters
    double T = std::atof(argv[2]), h = std::atof(argv[1]);
    unsigned int N = static_cast<unsigned int>(T / h);

    // Random seed
    std::default_random_engine generator{(unsigned int) time(NULL)};

    std::string outputFilename = std::string(DATA_PATH) + argv[3] + ".txt";
    std::ofstream output;
    output.open(outputFilename, std::ofstream::out | std::ofstream::trunc);

    // Reference solution
    RungeKutta refSolver(ODE, param, tableauRef);
    double hRef = h / 50;
    unsigned int NRef = static_cast<unsigned int>(T / hRef);
    std::vector<VectorXd> refSolution;
    VectorXd currSolution = ODE.initialCond;

    /* for (unsigned int i = 0; i < NRef; i++) {
        currSolution = refSolver.oneStep(currSolution, hRef);
        if ((i+1) % 50 == 0) {
            refSolution.push_back(currSolution);
        }
    }
    std::cout << refSolution.back().transpose() << std::endl; */

    // Initialization
    RungeKuttaRandomH Method(&generator, ODE,
                             param, tableauRef, h, p);

    VectorXd solution = ODE.initialCond;
    double error;

    double refHam = todaHamiltonian(ODE.initialCond);

    for (unsigned int i = 0; i < N; i++) {
        solution = Method.oneStep(solution);
        // error = (refSolution[i] - solution).norm();
        // output << h * (i+1) << "\t" << error << std::endl;
        output << std::fixed << std::setprecision(20) << h * (i+1) << "\t"
               << std::abs(todaHamiltonian(solution) - refHam) << std::endl;
    }

    output.close();

    return 0;
}
