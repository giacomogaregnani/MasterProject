#ifndef MODELINGERRORALL_HPP
#define MODELINGERRORALL_HPP

#include <EulerMaruyama.hpp>
#include <ParFilLib.hpp>

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
    std::vector<double> weights;
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
    void computePFAlt(unsigned int nMC);
    void getModErr(std::vector<std::vector<double>>& getData);
    void getWeights(std::vector<double>& getData);
};



#endif
