#include "Filter.hpp"

#include "../matplotlib-cpp-master/matplotlibcpp.h"
namespace plt = matplotlibcpp;


Filter::Filter(double delta, double beta, double h, unsigned long N) :
    N(N),
    h(h)
{
    filter.resize(N);
    Cb = beta / std::tgamma(1.0 / beta);

    for (unsigned long i = 0 ; i < N; i++) {
        // filter[i] = Cb / std::pow(delta, 1.0 / beta) * std::exp(-std::pow(h*i, beta) / delta);
        filter[i] = Cb / delta * std::exp(-std::pow(h*i / delta, beta));
    }
}

std::vector<double>& Filter::compute(std::vector<double>& data)
{
    filterData.resize(N);
    std::fill(filterData.begin(), filterData.end(), 0);

    for (unsigned int i = 0; i < N; i++) {
        for (unsigned int j = 0; j < i; j++) {
            filterData[i] += filter[i - j] * data[j];
        }
        filterData[i] *= h;
    }

    return filterData;
}

double Filter::getFilterConstant()
{
    return Cb;
}

double Filter::getFilterZero()
{
    return filter[0];
}