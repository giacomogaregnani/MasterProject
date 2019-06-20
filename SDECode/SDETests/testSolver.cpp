#include <EulerMaruyama.hpp>
#include <fstream>
#include <ctime>
#include <iostream>

double drift(double x, VectorXd& p)
{
    return -1.0 * p(0) * x;
}

double diffusion(double x, VectorXd &p)
{
    return std::sqrt(2.0 * p(1));
}

int main(int argc, char* argv[])
{
    oneDimSde sde;
    sde.drift = &drift;
    sde.diffusion = &diffusion;

    std::default_random_engine seed{2019};
    VectorXd param(2);
    param << 1.0, 0.5;

    EM1D solver(sde, param, seed);

    double T = 20;
    unsigned int N = 10000;
    auto h = T / N;
    unsigned int M = 1;
    double x, IC = 0;

    std::ofstream output(DATA_PATH + std::string("test.txt"), std::ofstream::out | std::ofstream::trunc);
    for (unsigned int i = 0; i < M; i++) {
        x = IC;
        output << 0.0 << "\t" << x << std::endl;
        for (unsigned int j = 0; j < N; j++) {
            x = solver.oneStep(h, x);
            output << (j+1)*h << "\t" << x << std::endl;
        }
    }
    output.close();

    return 0;
}