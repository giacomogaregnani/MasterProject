#ifndef FILTER_HPP
#define FILTER_HPP

#include <vector>
#include <cmath>

class Filter {
private:
    std::vector<double> filterData;
    std::vector<double> filter;
    unsigned long N;
    double h;
    double Cb;
public:
    Filter() = default;
    ~Filter() = default;
    Filter(double delta, double beta, double h, unsigned long N);
    std::vector<double>& compute(std::vector<double>& data);
    double getFilterConstant();
    double getFilterZero();
};


#endif //FILTER_HPP
