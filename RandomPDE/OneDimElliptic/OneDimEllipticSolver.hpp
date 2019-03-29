#ifndef ONEDIMELLSOLVER_H
#define ONEDIMELLSOLVER_H

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <memory>

using namespace Eigen;

class OneDimMesh {
private:
    // Extrema
    double xMin;
    double xMax;

    // Spacing
    double h;

    // Internal points
    VectorXd points;

    // Spacings
    VectorXd spacings;

public:
    OneDimMesh() {};

    OneDimMesh(double xMin, double xMax, double h);

    VectorXd& getPoints(void);

    void setPoints(VectorXd& pointsIn);

    VectorXd& getSpacings(void);
};

// Solve equation of type -div(kappa(x) grad(u)) = f(x);

class Solver {
private:

    // Coefficient
    double (*kappa) (double, VectorXd&);
    VectorXd kappaVec;
    bool isVectorKappa;

    // RHS
    double (*f) (double);

    // Mesh
    std::shared_ptr<OneDimMesh> mesh;
    size_t N;

    // Matrix
    SparseMatrix<double> A;

    // RHS-Vector
    VectorXd F;

    // Solution
    VectorXd U;

    // (Dirichlet) Boundary conditions
    double leftBC;
    double rightBC;

    // Assemblers
    void AssembleKappa(VectorXd& theta);
    void AssembleMatrix(void);
    void AssembleRHS(void);
    VectorXd Thomas();

public:
    Solver() {};

    Solver(double (*cond) (double, VectorXd&),
           double (*RHS) (double),
           double xMin, double xMax, double h,
           double leftBC, double rightBC);

    Solver(double (*RHS) (double),
           double xMin, double xMax, double h,
           double leftBC, double rightBC);

    VectorXd solve(VectorXd& theta);

    VectorXd evaluate(VectorXd& x);

    void changeMesh(VectorXd& points);

    VectorXd getMesh(void);
};


#endif // ONEDIMELLSOLVER_H