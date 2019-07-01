#ifndef MODELINGERRORALL_HPP
#define MODELINGERRORALL_HPP

#include <EulerMaruyama.hpp>
#include <ParFil.hpp>

class ModErrAll {
private:
    oneDimSde sdeCoarse;
    oneDimSde sdeFine;
    double (*V1) (double);
    double IC;
    VectorXd priorMean;
    VectorXd priorStdDev;
    double T;
    unsigned int N;
    std::vector<std::vector<double>> errors;
    std::vector<double> obs;
    double noise;

public:
    ModErrAll() = default;
    ~ModErrAll() = default;
    ModErrAll(oneDimSde sdeCoarse, oneDimSde sdeFine, double (*V1) (double),
           double IC, VectorXd &priorMean,
           VectorXd &priorStdDev, double T,
           unsigned int N, std::vector<double> &observations, double noise);
    void computePF(unsigned int nParam, unsigned int nMC);
    void getStats(std::vector<std::vector<double>>& getMean);
    VectorXd computeHomogeneous(VectorXd p, double L, double (*V1) (double));
};



#endif
