#include <Eigen/Dense>
#include <iostream>
#include <iomanip>

using namespace Eigen;

class RAMParam {
private:

    // RAM matrix
    MatrixXd S;

    // Desired acceptance ratio
    double targetAlpha;

    // Counter
    unsigned long count;

    // Size of the parameter space
    unsigned int nParam;

public:

    RAMParam() {};

    void init(double gamma, double alpha, unsigned int paramSize);

    void update(VectorXd& w, double alpha);

    MatrixXd& getS(void);
};


struct ProbDensFunc {
    double (*PDF) (VectorXd&);
};

struct LikelihoodFunc {
    double (*LIK) (VectorXd&, std::vector<VectorXd>&);
};


class MCMC {
private:

    // The proposal distribution
    std::normal_distribution<double> proposal;

    // The prior distribution
    ProbDensFunc* prior;

    // The likelihood function
    LikelihoodFunc* likelihood;

    // The probability generator
    std::uniform_real_distribution<double> probGen;

    // The number of samples
    unsigned long nMCMC;

    // The generated samples
    std::vector<VectorXd> samples;
    unsigned int sizeParam;

    // The type of update
    bool RAM;
    RAMParam RAMUpdate;

    // The observations
    std::vector<VectorXd> observations;

public:

    MCMC() {};

    MCMC(std::normal_distribution<double>::param_type& proposalParam,
         ProbDensFunc* prior, LikelihoodFunc* likelihood,
         std::vector<VectorXd>& observations,
         unsigned long nMCMC, bool RAM, VectorXd initGuess,
         double desiredAlpha);

    std::vector<VectorXd>& compute(std::default_random_engine* generator);
};