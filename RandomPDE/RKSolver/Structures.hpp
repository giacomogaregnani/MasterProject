#ifndef STRUCTURES_H
#define STRUCTURES_H

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>

using namespace Eigen;

// Collection of test problems;
enum problems {
    FITZNAG,
    LORENZ,
    HENHEIL,
    TEST1D,
    BRUSS,
    PEROX,
    KEPLER,
    KEPLERPERT,
    HIRES,
    PENDULUM,
    VDPOL,
    MODIFPEND
};

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
    GAUSS6,
    SYMPEULER,
    STORMVER
};

// Explicit or Implicit
enum ExplicitImplicit {
    IMPLICIT,
    EXPLICIT,
    STABEXP,
    SEPAR
};

// Definition of the ODE
struct odeDef {
    problems ode;
    int size;
    VectorXd initialCond;
    VectorXd (*odeFunc) (VectorXd, VectorXd&);
    MatrixXd (*odeJac) (VectorXd, VectorXd&);
    VectorXd (*exactSol) (VectorXd&, double);
    VectorXd refParam;
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