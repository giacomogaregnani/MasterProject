#include "UtilitiesGG.hpp"

VectorXd loadRefSolution(std::string& fileName, int solSize)
{
    VectorXd refSolution(solSize);

    std::fstream input;
    input.open(fileName, std::ofstream::in);
    for (int i = 0; i < solSize; i++) {
        input >> refSolution(i);
    }

    return refSolution;
}

void loadObservations(std::vector<double>& obsTimes, std::vector<VectorXd>& obs,
                      std::string& fileName, int solSize, double* noise)
{
    std::fstream input(fileName, std::ofstream::in);
    unsigned int nObs;
    input >> nObs;
    obsTimes.resize(nObs);
    obs.resize(nObs);

    for (unsigned int i = 0; i < nObs; i++) {
        input >> obsTimes[i];
        obs[i].resize(solSize);
        for (int j = 0; j < solSize; j++) {
            input >> obs[i](j);
        }
    }

    input >> *noise;
}

void computeOrderOfConvergence(std::vector<double>& errors,
                               std::vector<double>& orders,
                               double ratio)
{
    for (unsigned int j = 0; j < errors.size(); j++) {
        if (j > 0) {
            orders[j] = std::log(errors[j - 1] / errors[j]) / std::log(ratio);
        }
    }
}

void printConvergenceInfo(std::vector<double>& errors,
                          std::vector<double>& orders)
{
    double meanOrder = 0.0;

    std::cout << "========= CONVERGENCE INFO ===========" << std::endl;

    for (unsigned int j = 0; j < errors.size(); j++) {
        std::cout << "TIMESTEP " << j + 1 << " :\t";
        if (j > 0) {
            std::cout << "error = " << errors[j] << "\t";
            std::cout << "order = " << orders[j] << std::endl;
            meanOrder += orders[j];
        } else {
            std::cout << "error = " << errors[j] << std::endl;
        }
    }
    meanOrder /= orders.size();

    std::cout << "MEAN ORDER OF CONVERGENCE = " << meanOrder << std::endl;
    std::cout << "======================================" << std::endl;
}

VectorXd normalZeroMeanRandVec(int size,
                               std::default_random_engine* generator,
                               double stdDev)
{
    std::normal_distribution<double> normDist(0.0, stdDev);
    VectorXd randVec(size);

    for (int i = 0; i < size; i++) {
        randVec(i) = normDist(*generator);
    }

    return randVec;
}

std::vector<double> EigVecToStdVec(VectorXd& vec)
{
    std::vector<double> stdvec;
    for (unsigned int i = 0; i < vec.size(); i++) {
        stdvec.push_back(vec(i));
    }
    return stdvec;
}

VectorXd StdVecToEig(std::vector<double>& vec)
{
    VectorXd eigvec(vec.size());
    for (unsigned int i = 0; i < vec.size(); i++) {
        eigvec(i) = vec[i];
    }
    return eigvec;
}



