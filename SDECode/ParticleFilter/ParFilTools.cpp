#include "ParFilTools.hpp"
#include <iomanip>

double gaussianDensity(double& x, double& mu, double& sigma)
{
    return 1.0 / (std::sqrt(2.0 * M_PI) * sigma) * std::exp(-0.5 * (x  - mu) * (x - mu) / (sigma * sigma));
}

double gaussianDensity2d(Vector2d& x, Vector2d& mu, Matrix2d& sigma)
{
    double det = sigma.determinant();
    Vector2d diff = x - mu;
    LLT<Matrix2d> chol(sigma);
    return 1.0 / (2.0 * M_PI * std::sqrt(det)) * std::exp(-0.5 * diff.dot(chol.solve(diff)));
}

ForwardPFModErr::ForwardPFModErr(oneDimSde &sde, oneDimSde &sdeHomo, double h,
                                               std::default_random_engine &seed, double (*V1) (double)):
    sde(sde),
    sdeHomo(sdeHomo),
    seed(seed),
    h(h),
    V1(V1)
{
    solver = std::make_shared<EM1D>(sde, seed);
    solverHomo = std::make_shared<EM1D>(sdeHomo, seed);
    gaussian = std::normal_distribution<double>(0.0, std::sqrt(h));
    A(0, 0) = 1.0; A(0, 1) = 0.0; A(1, 0) = 1.0; A(1, 1) = -1.0;
}

std::vector<double> computeCoeffs(double (*V1) (double), double sigma, double L)
{
    unsigned int N = 10000;
    VectorXd discr;
    discr.setLinSpaced(N+1, 0.0, L);
    double h = L / N;

    std::vector<double> Zs = {0.0, 0.0};

    for (int i = 0; i < N; i++) {
        Zs[0] += h * (std::exp(V1(discr[i]) / sigma)  + std::exp(V1(discr[i+1])  / sigma)) / 2;
        Zs[1] += h * (std::exp(-V1(discr[i]) / sigma) + std::exp(-V1(discr[i+1]) / sigma)) / 2;
    }

    return Zs;
}

VectorXd ForwardPFModErr::computeHomogeneous(VectorXd param, double L, double (*V1) (double))
{
    param(1) = std::exp(param(1));
    param(2) = std::exp(param(2));

    VectorXd homParam(3);
    auto Zs = computeCoeffs(V1, param(2), L);

    homParam(1) = param(1) * L * L / (Zs[0] * Zs[1]);
    homParam(2) = param(2) * L * L / (Zs[0] * Zs[1]);

    homParam(1) = std::log(homParam(1));
    homParam(2) = std::log(homParam(2));

    homParam(0) = param(0);
    return homParam;
}

void ForwardPFModErr::modifyParam(VectorXd& theta)
{
    param = theta;
    VectorXd thetaHomo = computeHomogeneous(theta, 2*M_PI, V1);
    paramHomo = thetaHomo;
    solver->modifyParam(theta);
    solverHomo->modifyParam(thetaHomo);
}


Vector2d ForwardPFModErr::generateSample(Vector2d& oldValue)
{
    double BM = gaussian(seed);
    Vector2d newValue;
    newValue(0) = solver->oneStepGivenNoise(h, oldValue(0), BM);
    double temp = solverHomo->oneStepGivenNoise(h, oldValue(0) - oldValue(1), BM);
    newValue(1) = newValue(0) - temp;
    return newValue;
}

/* TODO: the IS approach tracks better the multiscale process and has a bigger ESS, but the moderr estimation looks
worse than the boostrap PF approach */

Vector2d ForwardPFModErr::generateSampleIS(Vector2d& oldValue, double newObs, double noise)
{
    Vector2d newValue;

    double alpha = sde.drift(oldValue(0), param),
           beta  = sde.diffusion(oldValue(0), param);

    // In Golightly, Wilkinson (2011) the diffusion is under square root
    beta *= beta;

    double denom = noise * noise + beta * h;
    double aj = alpha + beta * (newObs - (oldValue(0) + alpha * h)) / denom;
    double bj = beta - h * beta * beta / denom;

    ISmean(0) = oldValue(0) + aj * h;
    ISVariance(0, 0) = bj * h;

    double xHomOld = oldValue(0) - oldValue(1);
    ISmean(1) = xHomOld + h * sdeHomo.drift(xHomOld, paramHomo);
    double diffHom = sdeHomo.diffusion(xHomOld, paramHomo);
    ISVariance(1, 1) = diffHom * diffHom * h;
    ISVariance(1, 0) = 0.0;
    ISVariance(0, 1) = 0.0;

    gaussianIS.param(std::normal_distribution<double>::param_type(ISmean(0), std::sqrt(ISVariance(0, 0))));
    newValue(0) = gaussianIS(seed);
    gaussianIS.param(std::normal_distribution<double>::param_type(ISmean(1), std::sqrt(ISVariance(1, 1))));
    newValue(1) = newValue(0) - gaussianIS(seed);

    // We worked with x^\varepsilon and x^0 so we need to transform with A = (1 0; -1 1)
    ISmean = A * ISmean;
    ISVariance = A * ISVariance * A.transpose();

    return newValue;
}

double ForwardPFModErr::evalTransDensity(Vector2d &oldValue, Vector2d &newValue)
{
    double xHomOld = oldValue(0) - oldValue(1);
    transMean(0) = sde.drift(oldValue(0), param) * h + oldValue(0);
    transMean(1) = sdeHomo.drift(xHomOld , paramHomo) * h + xHomOld;

    double diffEps = sde.diffusion(oldValue(0), param),
           diffHom = sdeHomo.diffusion(xHomOld, paramHomo);
    double varEps = diffEps * diffEps * h,
           varHom = diffHom * diffHom * h,
           covEpsHom = diffEps * diffHom * h;

    transVariance(0, 0) = varEps;    transVariance(0, 1) = covEpsHom;
    transVariance(1, 0) = covEpsHom; transVariance(1, 1) = varHom;

    // We worked with x^\varepsilon and x^0 so we need to transform with A = (1 0; -1 1)
    transMean = A * transMean;
    transVariance = A * transVariance * A.transpose();

    // return gaussianDensity2d(newValue, transMean, transVariance); */
    double stddev = std::sqrt(transVariance(0, 0));
    return gaussianDensity(newValue(0), transMean(0), stddev);
}

double ForwardPFModErr::evalISDensity(Vector2d &newValue)
{
    double stddev = std::sqrt(ISVariance(0, 0));
    return gaussianDensity(newValue(0), ISmean(0), stddev);
}

