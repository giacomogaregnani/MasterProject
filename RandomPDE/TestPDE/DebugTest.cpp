#include <OneDimEllipticSolver.hpp>
#include <iostream>
#include <fstream>

#ifndef PI
#define PI 3.14159265359
#endif

double f(double x)
{
    return sin(2 * PI * x);
}

int main(int argc, char* argv[])
{
    VectorXd kappaVec = VectorXd::Ones(101);
    Solver solver(&f, 0.0, 1.0, 0.01, 0.0, 0.0);

    VectorXd theta;
    VectorXd u = solver.solve(kappaVec);

    std::fstream output(DATA_PATH + std::string("test.txt"), std::ofstream::out | std::ofstream::trunc);
    output << u.transpose() << std::endl;
    output.close();

    return 0;
}
