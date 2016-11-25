#include <vector>
#include <Eigen/Dense>
#include <Solver.hpp>

using namespace Eigen;

class Poisson {
private:
    static MatrixXd A;
public:
    static VectorXd initialCond;
    static VectorXd poissonTest(VectorXd argument, std::vector<double>& param);
    static MatrixXd poissonTestJ(VectorXd argument, std::vector<double>& param);
    static VectorXd poissonExact(std::vector<double>& param, double t);
    Poisson(int N);
};

VectorXd lorenz(VectorXd argument, std::vector<double>& param);
MatrixXd lorenzJ(VectorXd argument, std::vector<double>& param);
VectorXd testOneD(VectorXd argument, std::vector<double>& param);
MatrixXd testOneDJ(VectorXd argument, std::vector<double>& param);
VectorXd test(VectorXd argument, std::vector<double>& param);
VectorXd fitznag(VectorXd argument, std::vector<double>& param);
MatrixXd fitznagJ(VectorXd argument, std::vector<double>& param);
VectorXd vdPol(VectorXd argument, std::vector<double>& param);
MatrixXd vdPolJ(VectorXd argument, std::vector<double>& param);
VectorXd bruss(VectorXd argument, std::vector<double>& param);
MatrixXd brussJ(VectorXd argument, std::vector<double>& param);
VectorXd hires(VectorXd argument, std::vector<double>& param);
MatrixXd hiresJ(VectorXd argument, std::vector<double>& param);

void setProblem(odeDef* odeModel);
