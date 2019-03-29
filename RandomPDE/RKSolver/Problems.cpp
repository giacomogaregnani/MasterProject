#include <cmath>
#include <iostream>
#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>
#include "Tools.hpp"

#ifndef PI
#define PI 3.1415926535897
#endif

// PROBLEMS
VectorXd lorenz(VectorXd argument, VectorXd& param) {
    double sigma = param[0], rho = param[1], beta = param[2];

    VectorXd result(3);
    result(0) = sigma * (argument(1) - argument(0));
    result(1) = argument(0) * (rho - argument(2)) - argument(1);
    result(2) = argument(0) * argument(1) - beta * argument(2);
    return result;
}

MatrixXd lorenzJ(VectorXd argument, VectorXd& param) {
    double sigma = param[0], rho = param[1], beta = param[2];

    MatrixXd result(3, 3);

    result << - sigma, sigma, 0.0,
            rho - argument(2), - 1.0, - argument(0),
            argument(1), argument(0), - beta;

    return result;
}

VectorXd fitznag(VectorXd argument, VectorXd& param)
{
    double a = std::exp(param(0)),
           b = std::exp(param(1)),
           c = std::exp(param(2));
    VectorXd result(2);

    result(0) = c * (argument(0) - argument(0) * argument(0) * argument(0) / 3.0 + argument(1));
    result(1) = - 1.0 / c * (argument(0) - a + b * argument(1));

    return result;
}

MatrixXd fitznagJ(VectorXd argument, VectorXd& param)
{
    double b = param[1], c = param[2];
    MatrixXd result(2, 2);
    result << c * (1.0 - argument(0) * argument(0)), c,
            -1.0 / c, -1.0 * b / c;

    return result;
}

VectorXd HenHeil(VectorXd argument, VectorXd& param)
{
    VectorXd result(4);

    result << - argument(2) - 2.0 * argument(2) * argument(3),
            - argument(3) - argument(2) * argument(2)
            + argument(3) * argument(3),
            argument(0),
            argument(1);

    return result;
}

MatrixXd HenHeilJ(VectorXd argument, VectorXd& param)
{
    MatrixXd result = MatrixXd::Zero(4, 4);

    // p1
    result(0, 2) = - 1.0 - 2.0 * argument(3);
    result(0, 3) = - 2.0 * argument(2);

    // p2
    result(1, 2) = - 2.0 * argument(2);
    result(1, 3) = - 1.0 + 2.0 * argument(3);

    // q1
    result(2, 0) = 1.0;

    // q2
    result(3, 1) = 1.0;

    return result;
}

VectorXd test1D(VectorXd argument, VectorXd& param)
{
    return argument * param[0];
}

MatrixXd test1DJ(VectorXd argument, VectorXd& param)
{
    MatrixXd result(1, 1);
    result(0, 0) = param[0];

    return result;
}

VectorXd bruss(VectorXd argument, VectorXd& param)
{
    double alpha = param[0];
    int N = static_cast<int>(argument.size());
    int Nh = N / 2;
    double dX = 1.0 / (Nh + 1);
    double coeff = alpha / (dX * dX);
    VectorXd result(N);

    // u species
    result(0) = 1.0 + argument(0) * argument(0) * argument(Nh) - 4.0 * argument(0)
                + coeff * (1.0 - 2 * argument(0) + argument(1));
    for (int i = 1; i < Nh - 1; i++) {
        result(i) = 1.0 + argument(i) * argument(i) * argument(Nh + i)
                    - 4.0 * argument(i) + coeff * (argument(i - 1) - 2.0 * argument(i) + argument(i + 1));
    }
    result(Nh - 1) = 1.0 + argument(Nh - 1) * argument(Nh - 1) * argument(N - 1) - 4.0 * argument(Nh - 1)
                     + coeff * (argument(Nh - 2) - 2.0 * argument(Nh - 1) + 1.0);

    // v species
    result(Nh) = 3.0 * argument(0) - argument(0) * argument(0) * argument(Nh)
                 + coeff * (3.0 - 2.0 * argument(Nh) + argument(Nh + 1));
    for (int i = Nh + 1; i < N - 1; i++) {
        result(i) = 3.0 * argument(i - Nh) - argument(i - Nh) * argument(i - Nh) * argument(i)
                    + coeff * (argument(i - 1) - 2 * argument(i) + argument(i + 1));
    }
    result(N - 1) = 3.0 * argument(Nh - 1) - argument(Nh - 1) * argument(Nh - 1) * argument(N - 1)
                    + coeff * (argument(N - 2) - 2 * argument(N - 1) + 3.0);

    return result;
}

MatrixXd brussJ(VectorXd argument, VectorXd& param)
{
    double alpha = param[0];
    int N = static_cast<int>(argument.size());
    int Nh = N / 2;
    double dX = 1.0 / (Nh + 1);
    double coeff = alpha / (dX * dX);

    // Initialize the matrix
    MatrixXd result = MatrixXd::Zero(N, N);

    // du'_i / du_j and du'_i / dv_j
    result(0, 0) = 2.0 * argument(0) * argument(Nh) - 4.0 - 2.0 * coeff;
    result(0, 1) = coeff;
    result(0, Nh) = argument(0) * argument(0);
    for (int i = 1; i < Nh - 1; i++) {
        result(i, i) = 2.0 * argument(i) * argument(Nh + i) - 4.0 - 2.0 * coeff;
        result(i, i - 1) = coeff;
        result(i, i + 1) = coeff;
        result(i, i + Nh) = argument(i) * argument(i);
    }
    result(Nh - 1, Nh - 1) = 2.0 * argument(Nh - 1) * argument(N - 1) - 4.0 - 2.0 * coeff;
    result(Nh - 1, Nh - 2) = coeff;
    result(Nh - 1, N - 1) = argument(Nh - 1) * argument(Nh - 1);

    // dv'_i / du_j and dv'_i / dv_j
    result(Nh, 0) = 3.0 - 2.0 * argument(0) * argument(Nh);
    result(Nh, Nh) = -1.0 * argument(0) * argument(0) - 2.0 * coeff;
    result(Nh, Nh + 1) = coeff;
    for (int i = Nh + 1; i < N - 1; i++) {
        result(i, i) = -1.0 * argument(i - Nh) * argument(i - Nh) - 2.0 * coeff;
        result(i, i - 1) = coeff;
        result(i, i + 1) = coeff;
        result(i, i - Nh) = 3.0 - 2 * argument(i - Nh) * argument(i);
    }
    result(N - 1, N - 1) = -1.0 * argument(Nh - 1) * argument(Nh - 1) - 2 * coeff;
    result(N - 1, N - 2) = coeff;
    result(N - 1, Nh - 1) = 3.0 - 2 * argument(Nh - 1) * argument(N - 1);

    return result;
}

VectorXd PerOx(VectorXd argument, VectorXd& param)
{
    VectorXd result(4);

    double ABY = argument(0) * argument(1) * argument(3);
    double BX = argument(1) * argument(2);
    double XX = argument(2) * argument(2);
    result(0) = param[6] * (8.0 - argument(0)) - param[2] * ABY;
    result(1) = param[7] - param[0] * BX - param[2] * ABY;
    result(2) = param[0] * BX - 2.0 * param[1] * XX + 3 * param[2] * ABY - param[3] * argument(2) + param[5];
    result(3) = 2.0 * param[1] * XX - param[4] * argument(3) - param[2] * ABY;

    return result;
}

VectorXd Kepler(VectorXd argument, VectorXd& param)
{
    VectorXd result(4);

    double normQ = sqrt(argument(2) * argument(2) + argument(3) * argument(3));
    double normQCube = normQ * normQ * normQ;

    result << -argument(2) / normQCube,
            -argument(3) / normQCube,
            argument(0),
            argument(1);

    return result;
}

MatrixXd KeplerJ(VectorXd argument, VectorXd& param)
{
    MatrixXd result = MatrixXd::Zero(4, 4);

    double normQ = sqrt(argument(2) * argument(2) + argument(3) * argument(3));
    double normQSq = normQ * normQ;
    double normQV = normQSq * normQSq * normQ;

    // p1
    result(0, 2) = - (normQSq - 3.0 * argument(2) * argument(2)) / normQV;
    result(0, 3) = 3 * argument(2) * argument(3) / normQV;

    // p2
    result(1, 2) = result(0, 3);
    result(1, 3) = - (normQSq - 3.0 * argument(3) * argument(3)) / normQV;

    // q1
    result(2, 0) = 1.0;

    // q2
    result(3, 1) = 1.0;

    return result;
}

VectorXd KeplerPert(VectorXd argument, VectorXd& param)
{
    VectorXd result(4);

    double normQ = sqrt(argument(2) * argument(2) + argument(3) * argument(3));
    double normQSq = normQ * normQ;
    double normQCube = normQSq * normQ;
    double normQV = normQCube * normQSq;

    result << - argument(2) / normQCube - param[0] * argument(2) / normQV,
            -argument(3) / normQCube - param[0] * argument(3) / normQV,
            argument(0),
            argument(1);

    return result;
}

MatrixXd KeplerPertJ(VectorXd argument, VectorXd& param)
{
    MatrixXd result = MatrixXd::Zero(4, 4);

    double normQ = sqrt(argument(2) * argument(2) + argument(3) * argument(3));
    double normQSq = normQ * normQ;
    double normQV = normQSq * normQSq * normQ;

    // p1
    result(0, 2) = - (normQSq - 3.0 * argument(2) * argument(2)) / normQV
                   - param[0] * (normQSq - 5.0 * argument(2) * argument(2)) / (normQV * normQSq);
    result(0, 3) = 3 * argument(2) * argument(3) / normQV
                   + 5 * param[0] * argument(2) * argument(3) / (normQV * normQSq);

    // p2
    result(1, 2) = result(0, 3);
    result(1, 3) = - (normQSq - 3.0 * argument(3) * argument(3)) / normQV
                   - param[0] * (normQSq - 5.0 * argument(3) * argument(3)) / (normQV * normQSq);

    // q1
    result(2, 0) = 1.0;

    // q2
    result(3, 1) = 1.0;

    return result;
}

VectorXd Pendulum(VectorXd argument, VectorXd& param)
{
    VectorXd result(2);

    result(0) = -std::sin(argument(1));
    result(1) = argument(0);

    return result;
}

MatrixXd PendulumJ(VectorXd argument, VectorXd& param)
{
    MatrixXd result = MatrixXd::Zero(2,2);

    result(0, 1) = -std::cos(argument(1));
    result(1, 0) = 1.0;

    return result;
}

VectorXd modifPendulum(VectorXd argument, VectorXd& param)
{
    VectorXd result(2);

    result(0) = -std::sin(argument(1)) + 0.4 * std::cos(2.0 * argument(1));
    result(1) = argument(0);

    return result;
}

MatrixXd modifPendulumJ(VectorXd argument, VectorXd& param)
{
    MatrixXd result = MatrixXd::Zero(2,2);

    result(0, 1) = -std::cos(argument(1)) - 0.8 * std::sin(2.0 * argument(1));
    result(1, 0) = 1.0;

    return result;
}

VectorXd hires(VectorXd y, VectorXd& param)
{
    VectorXd result(8);

    result(0) = - param[0] * y(0) + 0.43 * y(1) + param[2] * y(2) + 0.0007;
    result(1) = param[0] * y(0) - 8.75 * y(1);
    result(2) = param[1] * y(2) + 0.43 * y(3) + 0.035 * y(4);
    result(3) = param[2] * y(1) + param[0] * y(2) - 1.12 * y(3);
    result(4) = param[3] * y(4) + 0.43 * y(5) + 0.43 * y(6);
    result(5) = - param[4] * y(5) * y(7) + 0.69 * y(3) + param[0] * y(4) -
                0.43 * y(5) + 0.69 * y(6);
    result(6) = param[4] * y(5) * y(7) - 1.81 * y(6);
    result(7) = - result(6);

    return result;
}

MatrixXd hiresJ(VectorXd y, VectorXd& param)
{
    MatrixXd result(8, 8);
    result << - param[0], 0.43, param[2], 0.0, 0.0, 0.0, 0.0, 0.0,
            param[0], - 8.75, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, param[1], 0.43, 0.035, 0.0, 0.0, 0.0,
            0.0, param[2], param[0], - 1.12, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, param[3], 0.43, 0.43, 0.0,
            0.0, 0.0, 0.0, 0.69, param[0], - 0.43 - param[4] * y(7), 0.69, - param[4] * y(5),
            0.0, 0.0, 0.0, 0.0, 0.0, param[4] * y(7), - 1.81, param[4] * y(5),
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1, 0.0;
    return result;
}

VectorXd vdPol(VectorXd y, VectorXd& param)
{
    VectorXd result(2);

    result(0) = y(1);
    result(1) = ((1 - y(0) * y(0)) * y(1) - y(0)) / param(0);

    return result;
}

MatrixXd vdPolJ(VectorXd y, VectorXd& param)
{
    MatrixXd result(2, 2);

    result(0, 0) = 0.0;
    result(0, 1) = 1.0;
    result(1, 0) = -1.0 / param(0) * (2.0 * y(0) * y(1) + 1.0);
    result(1, 1) = 1.0 / param(0) * (1.0 - y(0) * y(0));

    return result;
}


void setProblem(odeDef* odeModel)
{
    switch (odeModel->ode) {
        case HIRES:
            odeModel->size = 8;
            (odeModel->initialCond).resize(odeModel->size);
            for (int i = 1; i < 7; i++) {
                (odeModel->initialCond)(i) = 0.0;
            }
            (odeModel->initialCond)(0) = 1.0;
            (odeModel->initialCond)(7) = 0.0057;
            odeModel->odeFunc = &hires;
            odeModel->odeJac = &hiresJ;
            (odeModel->refParam).resize(5);
            odeModel->refParam(0) = 1.71;
            odeModel->refParam(1) = -10.03;
            odeModel->refParam(2) = 8.32;
            odeModel->refParam(3) = -1.745;
            odeModel->refParam(4) = 280.0;
            break;
        case FITZNAG:
            odeModel->size = 2;
            (odeModel->initialCond).resize(odeModel->size);
            (odeModel->initialCond)(0) = -1.0;
            (odeModel->initialCond)(1) = 1.0;
            odeModel->odeFunc = &fitznag;
            odeModel->odeJac = &fitznagJ;
            (odeModel->refParam).resize(3);
            odeModel->refParam(0) = std::log(0.2);
            odeModel->refParam(1) = std::log(0.2);
            odeModel->refParam(2) = std::log(3.0);
            break;
        case LORENZ:
            (odeModel->size) = 3;
            (odeModel->initialCond).resize(odeModel->size);
            (odeModel->initialCond)(0) = -10.0;
            (odeModel->initialCond)(1) = -1.0;
            (odeModel->initialCond)(2) = 40.0;
            odeModel->odeFunc = &lorenz;
            odeModel->odeJac = &lorenzJ;
            (odeModel->refParam).resize(3);
            (odeModel->refParam)(0) = 10.0;
            (odeModel->refParam)(1) = 28.0;
            (odeModel->refParam)(2) = 8.0 / 3.0;
            break;
        case HENHEIL:
            odeModel->size = 4;
            (odeModel->initialCond).resize(odeModel->size);
            (odeModel->initialCond)(0) = 0.0;
            (odeModel->initialCond)(1) = 0.1;
            (odeModel->initialCond)(2) = 0.5;
            (odeModel->initialCond)(3) = 0.0;
            odeModel->odeFunc = &HenHeil;
            odeModel->odeJac = &HenHeilJ;
            odeModel->refParam = VectorXd::Zero(0);
            break;
        case TEST1D:
            odeModel->size = 1;
            (odeModel->initialCond).resize(odeModel->size);
            (odeModel->initialCond)(0) = 1.0;
            odeModel->odeFunc = &test1D;
            odeModel->odeJac = &test1DJ;
            (odeModel->refParam).resize(1);
            (odeModel->refParam)(0) = -1.0;
            break;
        case PEROX:
            odeModel->size = 4;
            odeModel->initialCond.resize(odeModel->size);
            (odeModel->initialCond)(0) = 6.0;
            (odeModel->initialCond)(1) = 58.0;
            (odeModel->initialCond)(2) = 0.0;
            (odeModel->initialCond)(3) = 0.0;
            odeModel->odeFunc = &PerOx;
            odeModel->refParam.resize(8);
            (odeModel->refParam)[0] = 0.35;
            (odeModel->refParam)[1] = 250.0;
            (odeModel->refParam)[2] = 0.035;
            (odeModel->refParam)[3] = 20.0;
            (odeModel->refParam)[4] = 5.35;
            (odeModel->refParam)[5] = 1e-5;
            (odeModel->refParam)[6] = 0.1;
            (odeModel->refParam)[7] = 0.825;
            break;
        case BRUSS:
            odeModel->size = 60;
            (odeModel->initialCond).resize(odeModel->size);
            double x;
            for (int i = 0; i < odeModel->size / 2; i++) {
                x = static_cast<double>(i + 1) / (odeModel->size / 2.0 + 1);
                (odeModel->initialCond)(i) = 1.0 + sqrt(2.0) / (2.0 * PI) * sin(2.0 * PI * x);
            }
            for (int i = odeModel->size / 2; i < odeModel->size; i++) {
                (odeModel->initialCond)(i) = 3.0;
            }
            odeModel->odeFunc = &bruss;
            odeModel->odeJac = &brussJ;
            (odeModel->refParam).resize(1);
            (odeModel->refParam)[0] = 1.0 / 50.0;
            break;
        case KEPLER:
            odeModel->size = 4;
            (odeModel->initialCond).resize(odeModel->size);
            (odeModel->initialCond)(0) = 1.0 - 0.6;
            (odeModel->initialCond)(1) = 0.0;
            (odeModel->initialCond)(2) = 0.0;
            (odeModel->initialCond)(3) = std::sqrt((1.0 + 0.6) / (1.0 - 0.6));

            odeModel->odeFunc = &Kepler;
            odeModel->odeJac = &KeplerJ;
            odeModel->refParam = {};
            break;
        case KEPLERPERT:
            odeModel->size = 4;
            (odeModel->initialCond).resize(odeModel->size);
            (odeModel->initialCond)(0) = 0.0;
            (odeModel->initialCond)(1) = 2.0;
            (odeModel->initialCond)(2) = 0.4;
            (odeModel->initialCond)(3) = 0.0;

            odeModel->odeFunc = &KeplerPert;
            odeModel->odeJac = &KeplerPertJ;
            (odeModel->refParam).resize(1);
            (odeModel->refParam)[0] = 0.015;
            break;
        case PENDULUM:
            odeModel->size = 2;
            (odeModel->initialCond).resize(odeModel->size);
            (odeModel->initialCond)(0) = 1.5;
            (odeModel->initialCond)(1) = -PI;

            odeModel->odeFunc = &Pendulum;
            odeModel->odeJac = &PendulumJ;
            odeModel->refParam = {};
            break;
        case MODIFPEND:
            odeModel->size = 2;
            (odeModel->initialCond).resize(odeModel->size);
            (odeModel->initialCond)(0) = 1.5;
            (odeModel->initialCond)(1) = -PI;

            odeModel->odeFunc = &modifPendulum;
            odeModel->odeJac = &modifPendulumJ;
            odeModel->refParam = {};
            break;
        case VDPOL:
            odeModel->size = 2;
            (odeModel->initialCond).resize(odeModel->size);
            (odeModel->initialCond)(0) = 2.0;
            (odeModel->initialCond)(1) = 0.0;

            odeModel->odeFunc = &vdPol;
            odeModel->odeJac = &vdPolJ;
            odeModel->refParam.resize(1);
            odeModel->refParam[0] = 4;
            break;
    }
}