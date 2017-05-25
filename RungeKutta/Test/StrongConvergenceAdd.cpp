#include <RandomTimeStep.hpp>
#include <UtilitiesGG.hpp>

#include <iomanip>

int main(int argc, char* argv[])
{
    odeDef ODE;
    ODE.ode = FITZNAG;
    setProblem(&ODE);

    std::vector<double> param = {0.2, 0.2, 3.0};
    Butcher tableau(KU3, EXPLICIT, 0);

    double T = 1.0, h = 0.01;
    unsigned int nMC = 1000, nExp = 5;

    std::default_random_engine generator{(unsigned int) time(NULL)};

    // Output file
    std::string outputFilename = std::string(DATA_PATH) + argv[1] + ".txt";
    std::ofstream output;
    output.open(outputFilename, std::ofstream::out | std::ofstream::trunc);

    // Reference solution
    std::string refFileName = std::string(DATA_PATH) + argv[2] + ".txt";
    VectorXd refSolution = loadRefSolution(refFileName, ODE.size);

    // Errors and orders of convergence
    std::vector<double> errors(nExp);
    std::vector<double> orders(nExp - 1);

    // "Timestep" loop
    for (unsigned int j = 0; j < nExp; j++) {

        VectorXd MCMean = VectorXd::Zero(ODE.size);
        int N = static_cast<int>(T / h);

        std::cout << "Timestep = "
                  << h
                  << "\t\t"
                  << "Number of steps = "
                  << N
                  << std::endl;

        // MC loop
        unsigned int k;
        std::vector<double> MCError(nMC);

        #pragma omp parallel for num_threads(30) private(k)
        for (k = 0; k < nMC; k++) {

            RungeKuttaAddNoise Method(&generator, ODE, param, tableau, h, 3.0);
            VectorXd solution = ODE.initialCond;

            // Time integration loop
            for (int i = 0; i < N; i++) {
                solution = Method.oneStep(solution);
            }

            MCError[k] = (solution - refSolution).norm();
        }

        double strongError = 0.0;
        for (auto it : MCError) {
            strongError += it;
        }
        strongError /= nMC;

        errors[j] = strongError;
        output << std::fixed << std::setprecision(20) << h << "\t" << strongError << "\n";

        h /= 2;
    }

    computeOrderOfConvergence(errors, orders, 2.0);
    printConvergenceInfo(errors, orders);

    output.close();

    return 0;
}