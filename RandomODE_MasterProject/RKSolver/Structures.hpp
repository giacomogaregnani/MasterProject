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
    TEST,
    TEST1D,
    VDPOL,
    ROBERTSON,
    BRUSS,
    HEAT,
    HIRES,
    GAUSSMAP,
    LOGMAP,
    AUTOCAT,
    PEROX,
    KEPLER,
    KEPLERPERT,
    HENHEIL,
    TODA
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
    GAUSS4
};

// Explicit or Implicit
enum ExplicitImplicit {
    IMPLICIT,
    EXPLICIT,
    STABEXP
};

// Definition of the ODE
struct odeDef {
    problems ode;
    int size;
    VectorXd initialCond;
    VectorXd (*odeFunc) (VectorXd, std::vector<double>&);
    MatrixXd (*odeJac) (VectorXd, std::vector<double>&);
    VectorXd (*exactSol) (std::vector<double>&, double);
    std::vector<double> refParam;
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
};

void setProblem(odeDef* ODE, int nHeat = 0);

#endif //STRUCTURES_H