#include <MCMC.hpp>
#include <fstream>

double prior(VectorXd& x)
{
    double tmp = x(0) * x(0) - x(1);
    return -10 * tmp * tmp - pow(x(1) - 0.25, 4);
}

double likelihood(VectorXd& x, std::vector<VectorXd>& obs)
{
    return 0.0;
}

int main(int argc, char* argv[])
{
    if (argc < 2) {
        throw std::invalid_argument("Insert the desired accuracy\n");
    }

    // Init prior and likelihood
    ProbDensFunc priorObj;
    priorObj.PDF = &prior;
    LikelihoodFunc likObj;
    likObj.LIK = &likelihood;

    // Proposal parameter
    std::normal_distribution<double>::param_type param(0.0, 1.0);

    // Observations
    std::vector<VectorXd> obs = {};

    // MCMC Param
    unsigned long nMCMC = 100000;
    bool RAM = true;
    double desiredAlpha = std::atof(argv[1]);

    // First guess
    VectorXd initGuess(2);
    initGuess << 1.0, -1.0;

    // Initialize MCMC
    MCMC algorithm(param, &priorObj, &likObj, obs,
                   nMCMC, RAM, initGuess, desiredAlpha);
    std::default_random_engine generator{(unsigned int) time(NULL)};
    std::vector<VectorXd> result = algorithm.compute(&generator);


    // Output
    std::ofstream output(std::string(DATA_PATH) + "resultsMCMC.txt",
                         std::ofstream::out | std::ofstream::trunc);

    for (auto it : result) {
        output << it.transpose() << std::endl;
    }

    output.close();

    return 0;
}