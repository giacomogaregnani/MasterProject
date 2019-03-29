#ifndef PARTICLEFILTER_HPP
#define PARTICLEFILTER_HPP

#include <EulerMaruyama.hpp>
#include <random>
#include <memory>

// Boostrap particle filter

class ParFil {
private:
    std::shared_ptr<EM1D> Solver;
    double T;
    double IC;
    unsigned int samplingRatio;
    std::vector<double> obs;
    double noise;
    std::vector<std::vector<double>> X;
    unsigned long nParticles;
    double eps;
    std::vector<double> W;
    double likelihood;
    std::default_random_engine shufflerSeed;
    std::discrete_distribution<unsigned int> shuffler;

public:
    ParFil() = default;
    ~ParFil() = default;
    ParFil(std::vector<double>& y, double T, double IC,
           unsigned int sR, double noise, oneDimSde sde,
           double eps, unsigned long nParticles);
    void compute(VectorXd& theta);
    double getLikelihood() const;
    std::vector<double> getBestX();
};

#endif