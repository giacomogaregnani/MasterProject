#ifndef PROJECTTEST_STRUCTURES_H
#define PROJECTTEST_STRUCTURES_H

#include <Eigen/Dense>

using namespace Eigen;

enum stabMethods {
    stdRKC,
    ROCK
};

struct StabValues {
    double stiffIndex;
    double damping;
    int nStages;
    MatrixXd (*jacobian) (VectorXd, std::vector<double>&);
    stabMethods method;
};

enum problems {
    FITZNAG,
    LORENZ,
    TEST,
    TEST1D,
    VDPOL,
    ROBERTSON,
    BRUSS,
    POISSON,
    HIRES,
};

struct odeDef {
    problems ode;
    int size;
    VectorXd initialCond;
    VectorXd (*odeFunc) (VectorXd, std::vector<double>&);
    MatrixXd (*odeJac) (VectorXd, std::vector<double>&);
    VectorXd (*exactSol) (std::vector<double>&, double);
    std::vector<double> refParam;
};

enum implicitMethods {
    GAUSS,
    LOBATTO,
    RADAU
};

class Butcher {
private:
    implicitMethods method;
    MatrixXd A;
    VectorXd b;
    int stages;

public:
    Butcher(implicitMethods chosenMethod, int s);
    int getStages(void);
    MatrixXd getA(void);
    VectorXd getB(void);
};

#endif //PROJECTTEST_STRUCTURES_H
