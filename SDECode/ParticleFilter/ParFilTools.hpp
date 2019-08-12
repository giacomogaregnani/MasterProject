#ifndef PARFILTOOLS_HPP
#define PARFILTOOLS_HPP

#include <EulerMaruyama.hpp>
#include <memory>

double gaussianDensity(double& x, double& mu, double& sigma);
double gaussianDensity2d(Vector2d& x, Vector2d& mu, Matrix2d& sigma);

class ForwardPFModErr {
private:
    oneDimSde sde;
    oneDimSde sdeHomo;
    std::shared_ptr<EM1D> solver;
    std::shared_ptr<EM1D> solverHomo;
    double h;
    std::default_random_engine seed;
    double (*V1) (double);
    std::normal_distribution<double> gaussian;
    std::normal_distribution<double> gaussianIS;
    Vector2d ISmean;
    Matrix2d ISVariance;
    Vector2d transMean;
    Matrix2d transVariance;
    VectorXd param;
    VectorXd paramHomo;
    Matrix2d A;

public:
    ForwardPFModErr() = default;
    ~ForwardPFModErr() = default;
    ForwardPFModErr(oneDimSde& sde, oneDimSde& sdeHomo, double h,
                           std::default_random_engine& seed, double (*V1) (double));
    void modifyParam(VectorXd& theta);
    VectorXd computeHomogeneous(VectorXd param, double L, double (*V1) (double));
    Vector2d generateSample(Vector2d& oldValue);
    Vector2d generateSampleIS(Vector2d& oldValue, double newObs, double noise);
    double evalTransDensity(Vector2d& oldValue, Vector2d& newValue);
    double evalISDensity(Vector2d& newValue);
};

#endif