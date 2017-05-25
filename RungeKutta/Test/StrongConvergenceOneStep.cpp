#include <RandomTimeStep.hpp>
#include <UtilitiesGG.hpp>
#include <iomanip>

int main(int argc, char* argv[])
{
    odeDef ODE;
    ODE.ode = FITZNAG;
    setProblem(&ODE);
    Butcher tableau(RK4, EXPLICIT, 0);

    double h = 0.1, hRef = 1e-7;
    unsigned int nMC = 10000, nExp = 5;
    VectorXd refSolution(ODE.size);
    RungeKutta refRK(ODE, ODE.refParam, tableau);

    std::default_random_engine generator{(unsigned int) time(NULL)};

    std::vector<double> errors(nExp);
    std::vector<double> orders(nExp - 1);

    // "Timestep" loop
    for (unsigned int j = 0; j < nExp; j++) {

        std::cout << "Mean timestep = "
                  << h
                  << std::endl;

        // Reference solution
        std::cout << "Building ref solution ... " << std::endl;
        refSolution = ODE.initialCond;
        unsigned int NRef = static_cast<unsigned int>(h / hRef);
        for (unsigned int i = 0; i < NRef; i++) {
            refSolution = refRK.oneStep(refSolution, hRef);
        }

        // Monte Carlo loop
        std::cout << "Monte Carlo loop ... " << std::endl;
        VectorXd MCMean = VectorXd::Zero(ODE.size);
        unsigned int k;
        std::vector<double> finalTime(nMC), MCError(nMC);

        #pragma omp parallel for num_threads(30) private(k)
        for (k = 0; k < nMC; k++) {
            RungeKuttaRandomH Method(&generator, ODE, ODE.refParam, tableau, h, std::atof(argv[1]));
            VectorXd solution = Method.oneStep(ODE.initialCond);
            MCError[k] = (solution - refSolution).norm();
        }

        double strongError = 0.0;
        for (auto it : MCError) {
            strongError += it;
        }
        strongError /= nMC;
        errors[j] = strongError;

        // Time step update
        h /= 2;
    }

    // Error statistics
    computeOrderOfConvergence(errors, orders, 2);
    printConvergenceInfo(errors, orders);

    return 0;
}
