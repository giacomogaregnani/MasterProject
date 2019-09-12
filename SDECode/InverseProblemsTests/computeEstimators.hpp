#ifndef COMPUTEESTIMATORS_HPP
#define COMPUTEESTIMATORS_HPP

#include <vector>
#include <Eigen/Dense>

using namespace Eigen;

double estimateSigma(std::vector<double>&, double);

double estimateA(std::vector<double>&, double, double (*) (double));

double estimateA2(std::vector<double>&, std::vector<double>&, double, double (*) (double));

double estimateA3(std::vector<double>&, std::vector<double>&, double, double (*) (double));

std::vector<double> averageSequence(std::vector<double>& xIn, unsigned int windSize);

#endif