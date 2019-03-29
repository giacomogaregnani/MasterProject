#include <EulerMaruyama.hpp>
#include <fstream>
#include <ctime>
#include <iostream>

double drift(double x, VectorXd& p)
{
    return -1.0 * p(0) * (x * x * x - x);
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
    std::normal_distribution<double> BMGen(0.0, 1.0);

    double T = 20;

    unsigned int nBM = 10000;
    std::vector<double> BM(nBM+1);
    double hBM = T / nBM;
    BM[0] = 0.0;
    for (unsigned int i = 1; i < nBM+1; i++) {
        BM[i] = BM[i-1] + std::sqrt(hBM) * BMGen(seed);
    }

    VectorXd param(2);
    param << 1.0, 0.5;

    EM1D solver(sde, param, seed);


    std::ofstream output(DATA_PATH + std::string("testGN.txt"), std::ofstream::out | std::ofstream::trunc);

    for (unsigned int N = 100; N < nBM+1; N*=10) {
        auto h = T / N;
        auto ratio = static_cast<unsigned int>(std::round(nBM / N));
        double x = 0;
        output << 0.0 << "\t" << x << std::endl;
        for (unsigned int j = 0; j < N; j++) {
            x = solver.oneStepGivenNoise(h, x, BM[(j+1)*ratio] - BM[j*ratio]);
            output << (j+1)*h << "\t" << x << std::endl;
        }
    }

    output.close();

    return 0;
}