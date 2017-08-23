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

class MCMC {
private:

    // The proposal distribution
    std::normal_distribution<double> proposal;

    // The posterior distribution
    double (*posterior) (VectorXd& theta);

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


public:

    MCMC() {};

    MCMC(VectorXd& initGuess,
         std::normal_distribution<double>::param_type& proposalParam,
         double (*post) (VectorXd& theta),
         unsigned long nMCMC,
         bool RAM, double desiredAlpha);

    std::vector<VectorXd>& compute(std::default_random_engine* generator, bool noisy);
};