#include <MCMC.hpp>
#include <OneDimEllipticSolver.hpp>
#include <fstream>
#include <iomanip>

#ifndef PI
#define PI 3.14159265359
#endif

// Problem data

double f(double x)
{
    return 1.0;
    // return sin(2 * PI * x);
}

/* double kappaTheta(double x, VectorXd &theta)
{
    size_t nInt = static_cast<size_t> (theta.size());
    return exp(theta(static_cast<int>(x * nInt)));
} */

double kappaTheta(double x, VectorXd &theta)
{
    return theta(0);
}

// Posterior distribution

// Data
VectorXd observations;
VectorXd xObs;
double noise;
double h;
Solver inverseSolver;
Solver inverseProbSolver;
std::default_random_engine generator{(unsigned int) time(0)};

double prior(VectorXd& theta)
{
    return -0.5 * theta.dot(theta);
}

double posterior(VectorXd& theta)
{
    // Likelihood
    VectorXd uInv = inverseSolver.solve(theta);
    VectorXd dist = inverseSolver.evaluate(xObs) - observations;
    double likelihood = -0.5 / (noise * noise) * dist.dot(dist);

    double priorVal = prior(theta);

    return priorVal + likelihood;
}

double probPosterior(VectorXd& theta)
{
    size_t nMC = 100;
    std::vector<double> likelihoods(nMC);
    VectorXd mesh = inverseProbSolver.getMesh();
    long N = mesh.size() - 1;
    VectorXd intPoints = mesh.segment(1, N - 1);
    VectorXd newMesh;

    for (size_t k = 0; k < nMC; k++) {
        newMesh = mesh;
        newMesh.segment(1, N - 1) += h * h * 0.5 * VectorXd::Random(N - 1);
        inverseProbSolver.changeMesh(newMesh);
        VectorXd uInv = inverseProbSolver.solve(theta);
        inverseProbSolver.changeMesh(mesh);
        VectorXd dist = inverseProbSolver.evaluate(xObs) - observations;
        likelihoods[k] = -0.5 / (noise * noise) * dist.dot(dist);
    }
    inverseProbSolver.changeMesh(mesh);

    std::vector<double>::iterator maxLikIt = std::max_element(likelihoods.begin(), likelihoods.end());
    double maxLik = *maxLikIt;
    likelihoods.erase(maxLikIt, maxLikIt);

    double sum = 0;
    for (auto it : likelihoods) {
        sum += exp(it - maxLik);
    }
    double likelihood = maxLik + std::log(1 + sum) - std::log(static_cast<double>(nMC));
    double priorVal = prior(theta);

    return priorVal + likelihood;
}

double exactPosterior(VectorXd& theta)
{
    // Likelihood
    VectorXd exactSolution(xObs.size());
    for (int i = 0; i < xObs.size(); i++) {
        //exactSolution(i) = exp(-theta(0)) / (4*PI*PI) * sin(2 * PI * xObs(i));
        exactSolution(i) = 0.5 / theta(0) * (xObs(i) - xObs(i) * xObs(i));
    }
    VectorXd dist = exactSolution - observations;
    double likelihood = -0.5 / (noise * noise) * dist.dot(dist);

    double priorVal = prior(theta);

    return priorVal + likelihood;
}


int main(int argc, char* argv[])
{
    // INPUTS: h - nMCMC - noise stddev

    // Generate observations
    Solver refSolver(&kappaTheta, f, 0, 1, 1e-4, 0, 0);
    VectorXd theta(1);
    theta << 2;
    VectorXd uRef = refSolver.solve(theta);
    std::cout << theta.transpose() << std::endl;

    int nObs = 1;
    xObs.setLinSpaced(nObs, 0.25, 0.75);
    observations = refSolver.evaluate(xObs);

    noise = std::atof(argv[3]);
    std::normal_distribution<double> noiseDist(0.0, noise);

    for (int i = 0; i < nObs; i++) {
        observations(i) += noiseDist(generator);
    }

    // Set up the inverse solver
    h = std::atof(argv[1]);
    inverseSolver = Solver(&kappaTheta, f, 0, 1, h, 0, 0);
    inverseProbSolver = Solver(&kappaTheta, f, 0, 1, h, 0, 0);

    // MCMC
    VectorXd initGuess = VectorXd::Zero(theta.size());
    std::normal_distribution<double>::param_type proposalParam(0.0, 0.1);
    MCMC mcmc(initGuess, proposalParam, &posterior, std::stoul(argv[2]), true, 0.234);

    // Exact MCMC
    MCMC mcmcExact(initGuess, proposalParam, &exactPosterior, std::stoul(argv[2]), true, 0.234);

    // Prob MCMC
    MCMC mcmcProb(initGuess, proposalParam, &probPosterior, std::stoul(argv[2]), true, 0.234);

    // Solve the inverse problem
    std::vector<VectorXd> thetaAll = mcmc.compute(&generator, false);
    std::vector<VectorXd> thetaAllProb = mcmcProb.compute(&generator, false);
    std::vector<VectorXd> thetaAllExact = mcmcExact.compute(&generator, false);

    // Write on file
    std::fstream output(DATA_PATH + std::string("test.txt"), std::ofstream::out | std::ofstream::trunc);
    for (auto it : thetaAll) {
        output << std::fixed << std::setprecision(20) << it.transpose() << std::endl;
    }
    output.close();

    std::fstream outputExact(DATA_PATH + std::string("testExact.txt"), std::ofstream::out | std::ofstream::trunc);
    for (auto it : thetaAllExact) {
        outputExact << std::fixed << std::setprecision(20) << it.transpose() << std::endl;
    }
    outputExact.close();

    std::fstream outputProb(DATA_PATH + std::string("testProb.txt"), std::ofstream::out | std::ofstream::trunc);
    for (auto it : thetaAllProb) {
        outputProb << std::fixed << std::setprecision(20) << it.transpose() << std::endl;
    }
    outputProb.close();
    return 0;
}