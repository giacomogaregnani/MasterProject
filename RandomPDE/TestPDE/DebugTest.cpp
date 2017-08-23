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

double kappa(double x, VectorXd& theta)
{
    return 1.0;
}

int main(int argc, char* argv[])
{
    Solver solver(&kappa, &f, 0.0, 1.0, 0.1, 1.0, 2.0);

    VectorXd x = 0.5 * (VectorXd::Ones(10, 1) + VectorXd::Random(10, 1));
    std::sort(x.data(), x.data()+x.size());
    x(0) = 0; x(9) = 1;
    solver.changeMesh(x);

    VectorXd theta;
    VectorXd u = solver.solve(theta);

    std::fstream output(DATA_PATH + std::string("test.txt"), std::ofstream::out | std::ofstream::trunc);
    output << x.transpose() << std::endl << u.transpose() << std::endl;
    output.close();

    return 0;
}
