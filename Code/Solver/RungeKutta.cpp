#include <iostream>
#include <cmath>
#include "RungeKutta.hpp"

EulerForward::EulerForward(odeDef ODE, std::vector<double> paramVec)
{
	f = ODE.odeFunc;
	size = ODE.size;
	parameters = paramVec;
    theODE = ODE;
}

VectorXd EulerForward::oneStep(VectorXd solution, double h)
{
	return solution + h * f(solution, parameters);
}

int EulerForward::getOrder(void)
{
	return 1;
}

EulerBackwards::EulerBackwards(odeDef ODE, std::vector<double> paramVec)
{
    f = ODE.odeFunc;
    size = ODE.size;
	parameters = paramVec;
    theODE = ODE;
}

VectorXd EulerBackwards::oneStep(VectorXd solution, double h)
{
    return newton(theODE, solution, parameters, h, 100, 1e-12);
}

int EulerBackwards::getOrder(void)
{
	return 1;
}

RungeKutta::RungeKutta(odeDef ODE, std::vector<double> paramVec)
{
    f = ODE.odeFunc;
    size = ODE.size;
	parameters = paramVec;
    theODE = ODE;
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

MidPoint::MidPoint(odeDef ODE, std::vector<double> paramVec)
{
	f = ODE.odeFunc;
	size = ODE.size;
	parameters = paramVec;
	theODE = ODE;
}

VectorXd MidPoint::oneStep(VectorXd solution, double h)
{
	return solution + h * f(solution + 0.5 * h * f(solution, parameters), parameters);
}

int MidPoint::getOrder(void)
{
	return 2;
}

ImplicitRK::ImplicitRK(odeDef inODE, std::vector<double> paramVec,
					   MatrixXd inA, VectorXd inB, int inOrder)
{
	ODE = inODE;
	parameters = paramVec;
	A = inA;
	b = inB;
	order = inOrder;
}

VectorXd ImplicitRK::oneStep(VectorXd solution, double h)
{
	return newtonGeneral(ODE, solution, parameters, h, 100, 1e-12, A, b);
}

int ImplicitRK::getOrder(void)
{
	return order;
}