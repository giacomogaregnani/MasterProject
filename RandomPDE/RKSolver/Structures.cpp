#include "Structures.hpp"

Butcher::Butcher(methods chosenMethod, ExplicitImplicit type, int s)
{
    method = chosenMethod;
    methodType = type;
    switch (method) {
        case EULERFORWARD:
            stages = 1;
            A.resize(stages, stages);
            b.resize(stages);
            A(0, 0) = 0;
            b(0) = 1;
            break;
        case RK4:
            stages = 4;
            A.resize(stages, stages);
            b.resize(stages);
            for (int i = 0; i < stages; i++) {
                for (int j = i; j < stages; j++) {
                    A(i, j) = 0.0;
                }
            }
            A(1, 0) = 0.5;
            A(2, 0) = 0.0;
            A(2, 1) = 0.5;
            A(3, 0) = 0.0;
            A(3, 1) = 0.0;
            A(3, 2) = 1.0;

            b(0) = 1.0 / 6.0;
            b(1) = 1.0 / 3.0;
            b(2) = 1.0 / 3.0;
            b(3) = 1.0 / 6.0;
            break;
        case KU3:
            stages = 3;
            A.resize(stages, stages);
            b.resize(stages);
            for (int i = 0; i < stages; i++) {
                for (int j = i; j < stages; j++) {
                    A(i, j) = 0.0;
                }
            }
            A(1, 0) = 0.5;
            A(2, 0) = -1;
            A(2, 1) = 2.0;

            b(0) = 1.0 / 6.0;
            b(1) = 2.0 / 3.0;
            b(2) = 1.0 / 6.0;
            break;
        case EXPTRAPEZ:
            stages = 2;
            A.resize(stages, stages);
            b.resize(stages);
            A(0, 0) = 0.0;
            A(0, 1) = 0.0;
            A(1, 0) = 1.0;
            A(1, 1) = 0.0;

            b(0) = 0.5;
            b(1) = 0.5;
            break;
        case RKC:
            stages = s;
            break;
        case IMPMID:
            stages = 1;
            A.resize(stages, stages);
            b.resize(stages);
            A(0, 0) = 0.5;
            b(0) = 1.0;
            break;
        case IMPEULER:
            stages = 1;
            A.resize(stages, stages);
            b.resize(stages);
            A(0, 0) = 1.0;
            b(0) = 1.0;
            break;
        case IMPTRAPEZ:
            stages = 2;
            A.resize(stages, stages);
            b.resize(stages);
            A(0, 0) = 0.0;
            A(0, 1) = 0.0;
            A(1, 0) = 0.5;
            A(1, 1) = 0.5;
            b(0) = 0.5;
            b(1) = 0.5;
            break;
        case GAUSS4:
            stages = 2;

            A.resize(stages, stages);
            b.resize(stages);

            A(0, 0) = 0.25;
            A(0, 1) = 0.25 - sqrt(3) / 6.0;
            A(1, 0) = 0.25 + sqrt(3) / 6.0;
            A(1, 1) = 0.25;

            b(0) = 0.5;
            b(1) = 0.5;

            break;
        case GAUSS6:
            stages = 3;
            A.resize(stages, stages);
            b.resize(stages);
            b(0) = 5.0 / 18;
            b(1) = 4.0 / 9;
            b(2) = 5.0 / 18;
            A(0, 0) = 5.0 / 36;
            A(0, 1) = 2.0 / 9 - sqrt(15.0) / 15;
            A(0, 2) = 5.0 / 36 - sqrt(15.0) / 30;
            A(1, 0) = 5.0 / 36 + sqrt(15.0) / 24;
            A(1, 1) = 2.0 / 9;
            A(1, 2) = 5.0 / 36 - sqrt(15.0) / 24;
            A(2, 0) = 5.0 / 36 + sqrt(15.0) / 30;
            A(2, 1) = 2.0 / 9 + sqrt(15.0) / 15;
            A(2, 2) = 5.0 / 36;

            break;
        case SYMPEULER:
            break;
        case SSTAGETRAPEZ:
            stages = 5;
            A.resize(stages, stages);
            b.resize(stages);
            b(0) = 7.0;
            b(1) = 32.0;
            b(2) = 12.0;
            b(3) = 32.0;
            b(4) = 7.0;
            b /= 90.0;
            std::vector<double> coeff = {0.25, 0.5, 0.75, 1.0};

            for (int i = 0; i < stages; i++) {
                A(0, i) = 0.0;
            }

            for (int i = 1; i < stages; i++) {
                for (int j = 0; j < stages; j++) {
                    A(i, j) = b(j) * coeff[i - 1];
                }
            }
            break;
    }
}

MatrixXd Butcher::getA(void)
{
    return A;
}

VectorXd Butcher::getB(void)
{
    return b;
}

int Butcher::getStages(void)
{
    return stages;
}

ExplicitImplicit Butcher::getType(void)
{
    return methodType;
}

methods Butcher::getMethod(void)
{
    return method;
}