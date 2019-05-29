#ifndef MODELINGERROR_HPP
#define MODELINGERROR_HPP

#include <EulerMaruyama.hpp>
#include <ParFil.hpp>

class ModErr {
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
    std::vector<double> means;
    std::vector<double> stdDevs;
    std::vector<double> obs;
    double noise;

public:
    ModErr() = default;
    ~ModErr() = default;
    ModErr(oneDimSde sdeCoarse, oneDimSde sdeFine, double (*V1) (double),
           double IC, VectorXd &priorMean,
           VectorXd &priorStdDev, double T,
           unsigned int N, std::vector<double> &observations, double noise);
    void computePF(unsigned int nParam, unsigned int nMC);
    void compute(unsigned int nParam, unsigned int nMC);
    void getStats(std::vector<double>& getMean, std::vector<double>& getStdDev);
    VectorXd computeHomogeneous(VectorXd p, double L, double (*V1) (double));
};



#endif
