#include <iostream>
#include <cmath>
#include "RungeKutta.hpp"

EulerForward::EulerForward(int n, VectorXd (*func) (VectorXd, std::vector<double>&), std::vector<double> paramVec)
{
	f = func;
	size = n;
	parameters = paramVec;
}

VectorXd EulerForward::oneStep(VectorXd solution, double h)
{
	return solution + h * f(solution, parameters);                        
}

int EulerForward::getOrder(void)
{
	return 1;
}

RungeKutta::RungeKutta(int n, VectorXd (*func) (VectorXd, std::vector<double>&), std::vector<double> paramVec)
{
	f = func;
	size = n;
	parameters = paramVec;
}

VectorXd RungeKutta::oneStep(VectorXd solution, double h)
{
	VectorXd k1(size);
	VectorXd k2(size);
	VectorXd k3(size);
	VectorXd k4(size);

	k1 = f(solution, parameters);
	k2 = f(solution + h * 0.5 * k1, parameters);
	k3 = f(solution + h * 0.5 * k2, parameters);
	k4 = f(solution + h * k3, parameters);

	return solution + h * (1.0 / 6.0 * k1
                            + 1.0 / 3.0 * (k2 + k3)
                            + 1.0 / 6.0 * k4);
}

int RungeKutta::getOrder(void)
{
	return 4;
}
