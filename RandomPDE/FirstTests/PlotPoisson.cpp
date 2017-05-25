#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <stdexcept>
#include <Eigen/Sparse>
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;
typedef Triplet<double> T;

#ifndef PI
#define PI 3.14159265358979323846
#endif

VectorXd rightHandSide(double (*f) (double), VectorXd x)
{
    VectorXd F(x.size());
    double midPointOne, midPointTwo;
    long s = x.size();

    F(0) = x(0) * f(x(0) / 2) + (x(1) - x(0)) * f((x(1) + x(0)) / 2);

    for (int i = 1; i < x.size() - 1; i++) {
        midPointOne = x(i) + x(i-1);
        midPointTwo = x(i+1) + x(i);
        F(i) = (x(i) - x(i-1)) * f(midPointOne / 2)
               + (x(i+1) - x(i)) * f(midPointTwo / 2);
    }

    F(s-1) = (x(s) - x(s-1)) * f((x(s) + x(s-1)) / 2)
             + (1 - x(s)) * f((1.0 + x(s)) / 2);

    return F * 0.5;
}

double sinPiX(double x)
{
    return sin(PI * x);
}

double cosPiX(double x)
{
    return cos(PI * x);
}

int main (int argc, char* argv[])
{
    if (argc < 4) {
        throw invalid_argument("\nProgram arguments : \n"
                                       "1 : mesh size\n"
                                       "2 : number of MC\n"
                                       "3 : exponent p");
    }

    double h = atof(argv[1]), p = atof(argv[3]);
    unsigned long nMC = std::stoul(argv[2]);

    double (*f)(double) = &sinPiX;

    ofstream output(DATA_PATH + string("plotFEM.txt"), ofstream::out | ofstream::trunc);
    output << h << endl;

    int N = static_cast<int> (1.0 / h) - 1;

    SparseMatrix<double> A(N, N);
    vector<T> triplets;
    triplets.push_back(T(0, 0, 2.0));
    triplets.push_back(T(0, 1, -1.0));
    for (int i = 1; i < N - 1; i++) {
        triplets.push_back(T(i, i, 2.0));
        triplets.push_back(T(i, i-1, -1.0));
        triplets.push_back(T(i, i+1, -1.0));
    }
    triplets.push_back(T(N-1, N-1, 2.0));
    triplets.push_back(T(N-1, N-2, -1.0));
    A.setFromTriplets(triplets.begin(), triplets.end());
    A /= h;
    SimplicialLLT<SparseMatrix<double>> Cholesky(A);

    VectorXd x;
    VectorXd ones = VectorXd::Ones(N);
    x.setLinSpaced(N, h, 1.0 - h);

    vector<VectorXd> solutions(nMC);

    unsigned long i;
    #pragma omp parallel for private(i) num_threads(30)
    for (i = 0; i < nMC; i++) {
        VectorXd H(N), X(N);

        H.setRandom();
        X = x + H * pow(h, p);
        VectorXd F = rightHandSide(f, X);

        solutions[i] = Cholesky.solve(F);
    }

    VectorXd meanSolution(N);
    for (auto it : solutions) {
        output << 0 << "\t" << it.transpose() << "\t" << 0 << endl;
        meanSolution += it;
    }
    meanSolution /= nMC;
    output << 0 << "\t" << meanSolution.transpose() << "\t" << 0 << endl;

    output.close();

    return 0;
}