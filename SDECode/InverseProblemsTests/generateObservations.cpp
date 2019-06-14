#include "generateObservations.hpp"

std::vector<double> generateObservations1D(oneDimSde sde, double IC, VectorXd &param,
                                           double T, unsigned int N, unsigned int seedN)
{
    // std::random_device dev;
    std::default_random_engine seed{seedN};
    EM1D solver(sde, param, seed);

    auto h = T / N;
    std::vector<double> solution(N+1);
    solution[0] = IC;

    for (unsigned int j = 0; j < N; j++) {
        solution[j+1] = solver.oneStep(h, solution[j]);
    }

    return solution;
}