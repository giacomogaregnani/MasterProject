#ifndef STRUCTURES_H
#define STRUCTURES_H

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>

using namespace Eigen;

// Collection of classic RK methods
enum methods {
    EULERFORWARD,
    EXPTRAPEZ,
    KU3,
    RK4,
    RKC,
    IMPMID,
    IMPEULER,
    IMPTRAPEZ,
    SSTAGETRAPEZ,
    GAUSS4,
    GAUSS6
};

// Explicit or Implicit
enum ExplicitImplicit {
    IMPLICIT,
    EXPLICIT,
    STABEXP
};

class LinearEquation {
private:
    MatrixXd A;
    MatrixXd M;
    int size;
    VectorXd initialCond;
    FullPivLU<MatrixXd> LUMat;

public:
    VectorXd f(VectorXd);
    int getSize() {return size;};
    MatrixXd& getA() {return A};
    MatrixXd& getM() {return M};
};

// Butcher table class
class Butcher {
private:
    methods method;
    ExplicitImplicit methodType;
    MatrixXd A;
    VectorXd b;
    int stages;

public:
    Butcher() {};
    Butcher(methods chosenMethod, ExplicitImplicit type, int s = 0);
    int getStages(void);
    MatrixXd getA(void);
    VectorXd getB(void);
    ExplicitImplicit getType(void);
    methods getMethod(void);
};

void setProblem(odeDef* ODE);

#endif //STRUCTURES_H