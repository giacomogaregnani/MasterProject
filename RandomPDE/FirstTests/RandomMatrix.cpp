#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <stdexcept>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <Eigen/Dense>
#include <random>

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

SparseMatrix<double> assembleMatrix(VectorXd x)
{
    int N = x.size() - 1;
    SparseMatrix<double> A(N-1, N-1);

    A.insert(0, 0) = 1.0 / (x(1) - x(0)) + 1.0 / (x(2) - x(1));
    A.insert(0, 1) = -1.0 / (x(2) - x(1));
    for (int i = 1; i < N-2; i++) {
        A.insert(i, i) = 1.0 / (x(i+1) - x(i)) + 1.0 / (x(i+2) - x(i+1));
        A.insert(i, i+1) = -1.0 / (x(i+2) - x(i+1));
        A.insert(i, i-1) = -1.0 / (x(i+1) - x(i));
    }
    int i = N-2;
    A.insert(i, i) = 1.0 / (x(i+1) - x(i)) + 1.0 / (x(i+2) - x(i+1));
    A.insert(i, i-1) = -1.0 / (x(i+1) - x(i));

    return A;
}

double sinPiX(double x)
{
    return sin(PI * x);
}

double sinPiXUPrimeEx(double x)
{
    return cos(PI * x) / PI;
}

double poly(double x)
{
    return x;
}

double polyPrimeEx(double x)
{
    return - x * x / 2.0 + 1.0 / 6.0;
}

double FEMFirstDerivative(double yLeft, double yRight, double xLeft, double xRight)
{
    return (yLeft - yRight) / (xLeft - xRight);
}


double H1Error(double (*uPrimeEx) (double), VectorXd& FEMSol, VectorXd x, double h)
{
    unsigned int nQuad = 10;
    double integralTmp, integralTot = 0.0,
            hQuad = h / nQuad, femValue, tmp;

    // Inside elements
    for (int i = 1; i < x.size(); i++) {
        integralTmp = 0.0;
        femValue = FEMFirstDerivative(FEMSol(i), FEMSol(i-1), x(i), x(i-1));
        for (unsigned int j = 0; j < nQuad; j++) {
            tmp = femValue - uPrimeEx(x(i) + j * hQuad);
            integralTmp += tmp * tmp;
        }
        integralTot += integralTmp * hQuad;
    }

    return sqrt(integralTot);
}

int main (int argc, char* argv[])
{
    if (argc < 5) {
        throw invalid_argument("\nProgram arguments : \n"
                                       "1 : mesh size\n"
                                       "2 : number of MC\n"
                                       "3 : number of time step reductions\n"
                                       "4 : exponent p");
    }

    double h = atof(argv[1]), p = atof(argv[4]);
    unsigned long nMC = std::stoul(argv[2]);
    int nEx = atoi(argv[3]);

    vector<double> strongErrors((unsigned int) nEx, 0.0);
    vector<double> weakErrors((unsigned int) nEx, 0.0);

    double (*f)(double) = &poly;
    double (*uExPrime)(double) = &polyPrimeEx;

    default_random_engine generator{(unsigned int) time(NULL)};

    ofstream output(DATA_PATH + string("convFEM.txt"), ofstream::out | ofstream::trunc);

    vector<VectorXd> solutions(nMC);

    for (int k = 0; k < nEx; k++) {

        uniform_real_distribution<double> uniform(-pow(h, p), pow(h, p));

        // Number of elements
        int N = static_cast<int> (1.0 / h);

        // Create "base" mesh
        VectorXd x;
        x.setLinSpaced(N+1, 0.0, 1.0);

        unsigned long i;
        #pragma omp parallel for private(i) num_threads(30)
        for (i = 0; i < nMC; i++) {

            // Create random mesh
            VectorXd H(N+1), X(N+1);
            for (int i = 0; i < N+1; i++) {
                X(i) = x(i) + uniform(generator);
            }

            // Assemble matrix and r.h.s.
            VectorXd F = rightHandSide(f, X.segment(1, N-1));
            SparseMatrix<double> A = assembleMatrix(X);
            SimplicialLLT<SparseMatrix<double>> LU(A);

            // Solve
            solutions[i] = VectorXd::Zero(N+1);
            solutions[i].segment(1, N-1) = LU.solve(F);

            // Update strong error
            strongErrors[k] += H1Error(uExPrime, solutions[i], x, h);
        }
        strongErrors[k] /= nMC;

        VectorXd averageSolution = VectorXd::Zero(N+1);
        for (auto it : solutions) {
            averageSolution += it;
        }
        averageSolution /= nMC;
        weakErrors[k] = H1Error(uExPrime, averageSolution, x, h);

        cout << fixed << setprecision(6)
             << "h = " << h << " || strong H1 error = " << strongErrors[k]
             << " || weak H1 error " << weakErrors[k] << endl;

        output << h << "\t"
               << fixed << setprecision(18)
               << strongErrors[k] << "\t" << weakErrors[k] << endl;

        h /= 2;
    }

    output.close();

    return 0;
}