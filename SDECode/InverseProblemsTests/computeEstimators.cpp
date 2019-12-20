#include <iostream>
#include "computeEstimators.hpp"

double estimateSigma(std::vector<double>& x, double del)
{
    double sigma = 0.0, diff;
    auto N = x.size()-1;

    for (unsigned long i = 0; i < N; i++) {
        diff = x[i+1] - x[i];
        sigma += diff * diff;
    }

    sigma /= (2.0 * del * N);
    return sigma;
}

double estimateATilde(std::vector<double>& x, double del, double (*gradV) (double), double (*laplV) (double), double Sigma)
{
    double num = 0.0, denom = 0.0, temp;
    auto N = x.size() - 1;

    for (unsigned long i = 0; i < N; i++) {
        temp = gradV(x[i]);
        denom += temp * temp;
        num += laplV(x[i]);
    }

    return Sigma * del * num / (del * denom);
}

double estimateA(std::vector<double>& x, double del, double (*gradV) (double))
{
    double num = 0.0, denom = 0.0, temp;
    auto N = x.size() - 1;

    for (unsigned long i = 0; i < N; i++) {
        temp = gradV(x[i]);
        num += temp * (x[i+1] - x[i]);
        denom += temp * temp;
    }

   return -1.0 * num / (del * denom);
}

std::vector<double> averageSequence(std::vector<double>& xIn, unsigned int windSize)
{
    std::vector<double> xAvg(xIn.size());
    xAvg[0] = xIn[0];

    for (unsigned int i = 1; i < windSize; i++) // xAvg[i] = (i*xAvg[i-1] + xIn[i]) / (i+1);
        xAvg[i] = xAvg[i-1] + (xIn[i] - xIn[0]) / windSize;
    for (unsigned int i = windSize; i < xIn.size(); i++)
        xAvg[i] = xAvg[i-1] + (xIn[i] - xIn[i-windSize]) / windSize;

    return xAvg;
}
