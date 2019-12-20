#ifndef COMPUTEESTIMATORS_HPP
#define COMPUTEESTIMATORS_HPP

#include <vector>
#include <Eigen/Dense>

using namespace Eigen;

double estimateSigma(std::vector<double>&, double);

double estimateATilde(std::vector<double>& x, double del, double (*gradV) (double), double (*laplV) (double), double Sigma);

double estimateA(std::vector<double>&, double, double (*) (double));

std::vector<double> averageSequence(std::vector<double>& xIn, unsigned int windSize);

#endif