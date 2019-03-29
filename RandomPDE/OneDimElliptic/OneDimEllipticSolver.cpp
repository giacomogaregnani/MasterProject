#include <iostream>
#include <vector>
#include "OneDimEllipticSolver.hpp"

Solver::Solver(double (*kappa) (double, VectorXd&),
               double (*f) (double),
               double xMin, double xMax, double h,
               double leftBC, double rightBC):
        kappa(kappa),
        isVectorKappa(false),
        f(f),
        leftBC(leftBC),
        rightBC(rightBC)
{
    mesh = std::make_shared<OneDimMesh>(xMin, xMax, h);
    N = static_cast<size_t>(mesh->getPoints().size()) - 1;
}

Solver::Solver(double (*f) (double),
               double xMin, double xMax, double h,
               double leftBC, double rightBC):
        isVectorKappa(true),
        f(f),
        leftBC(leftBC),
        rightBC(rightBC)
{
    mesh = std::make_shared<OneDimMesh>(xMin, xMax, h);
    N = static_cast<size_t>(mesh->getPoints().size()) - 1;
}

void Solver::AssembleKappa(VectorXd& theta)
{
    if (isVectorKappa) {
        if (theta.size() == 1) {
            kappaVec = theta(0) * VectorXd::Ones(N+1);
        } else {
            kappaVec = theta;
        }
    }
    else {
        kappaVec.resize(N+1);

        VectorXd X = mesh->getPoints();

        for (size_t i = 0; i < N+1; i++) {
            kappaVec(i) = kappa(X(i), theta);
        }
    }
}

void Solver::AssembleMatrix(void)
{
    A.resize(N-1, N-1);
    VectorXd H = mesh->getSpacings();
    std::vector<Triplet<double>> tripletVector;
    double avg, avgPlusOne;

    avg = (kappaVec(0) + kappaVec(1)) / 2;
    avgPlusOne = (kappaVec(1) + kappaVec(2)) / 2;
    tripletVector.push_back(Triplet<double>(0, 0, 1.0 / H(0) * avg +
                                                  1.0 / H(1) * avgPlusOne));
    tripletVector.push_back(Triplet<double>(0, 1, -1.0 / H(1) * avgPlusOne));

    for (int i = 1; i < static_cast<int>(N)-2; i++) {
        tripletVector.push_back(Triplet<double>(i, i-1, tripletVector.back().value()));
        avg = avgPlusOne;
        avgPlusOne = (kappaVec(i) + kappaVec(i+1)) / 2;
        tripletVector.push_back(Triplet<double>(i, i, 1.0 / H(i) * avg +
                                                      1.0 / H(i+1) * avgPlusOne));
        tripletVector.push_back(Triplet<double>(i, i+1, -1.0 / H(i+1) * avgPlusOne));
    }

    int i = static_cast<int>(N) - 2;
    avg = avgPlusOne;
    avgPlusOne = (kappaVec(i) + kappaVec(i+1)) / 2;
    tripletVector.push_back(Triplet<double>(N-2, N-3, tripletVector.back().value()));
    tripletVector.push_back(Triplet<double>(N-2, N-2, 1.0 / H(N-2) * avg +
                                                      1.0 / H(N-1) * avgPlusOne));
    A.setFromTriplets(tripletVector.begin(), tripletVector.end());
}

void Solver::AssembleRHS(void)
{
    F.resize(N-1);

    VectorXd H = mesh->getSpacings();
    VectorXd X = mesh->getPoints();

    F(0) = f(X(1)) * (H(0) + H(1)) / 2.0 + leftBC * kappaVec(0) / H(0);
    F(N-2) = f(X(N-1)) * (H(N-2) + H(N-1)) / 2.0 + rightBC * kappaVec(N-1) / H(N-1);

    for (size_t i = 1; i < N-2; i++) {
        F(i) = f(X(i+1)) * (H(i) + H(i+1)) / 2.0;
    }
}

VectorXd Solver::Thomas()
{
    VectorXd c(N-2), d(N-1);
    double denom;

    c(0) = A.coeff(0, 1) / A.coeff(0, 0);
    d(0) = F(0) / A.coeff(0, 0);

    for (size_t i = 1; i < N-2; i++) {
        denom = A.coeff(i, i) - A.coeff(i, i-1) * c(i-1);
        c(i) = A.coeff(i, i+1) / denom;
        d(i) = (F(i) - A.coeff(i, i-1) * d(i-1)) / denom;
    }
    denom = A.coeff(N-2, N-2) - A.coeff(N-2, N-3) * c(N-3);
    d(N-2) = (F(N-2) - A.coeff(N-2, N-3) * d(N-3)) / denom;

    VectorXd u(N-1);

    u(N-2) = d(N-2);
    for (int i = N-3; i > -1; i--) {
        u(i) = d(i) - c(i) * u(i+1);
    }

    return u;
}

VectorXd Solver::solve(VectorXd& theta)
{
    AssembleKappa(theta);
    AssembleMatrix();
    AssembleRHS();
    U.resize(N+1);

    // Solve on internal nodes
    U.segment(1, N-1) = Thomas();

    // Apply BC
    U(0) = leftBC;
    U(N) = rightBC;

    return U;
}

VectorXd Solver::evaluate(VectorXd& x)
{
    long nObs = x.size();
    VectorXd values(nObs);
    VectorXd X = mesh->getPoints();
    VectorXd H = mesh->getSpacings();

    double span = X(N) - X(0);
    size_t n;

    for (long i = 0; i < nObs; i++) {
        n = static_cast<size_t>(x(i) * N / span);
        values(i) = U(n) + (U(n+1) - U(n)) / H(n) * (x(i) - X(n));
    }

    return values;
}

void Solver::changeMesh(VectorXd &points)
{
    mesh->setPoints(points);
    N = static_cast<size_t>(mesh->getPoints().size()) - 1;
}

VectorXd Solver::getMesh()
{
    return mesh->getPoints();
}