#include <vector>
#include <Eigen/Dense>

using namespace Eigen;

enum problems {
	FITZNAG,
	LORENZ,
	TEST,
	TEST1D,
	VDPOL,
	ROBERTSON,
	BRUSS,
	POISSON
};

class Poisson {
private:
    static MatrixXd A;
public:
    static VectorXd poissonTest(VectorXd argument, std::vector<double>& param);
    Poisson(int N);
};

VectorXd lorenz(VectorXd argument, std::vector<double>& param);
VectorXd testOneD(VectorXd argument, std::vector<double>& param);
VectorXd test(VectorXd argument, std::vector<double>& param);
VectorXd fitznag(VectorXd argument, std::vector<double>& param);
VectorXd vdPol(VectorXd argument, std::vector<double>& param);
VectorXd bruss(VectorXd argument, std::vector<double>& param);

void setProblem(VectorXd* initialCond, VectorXd (**odeFunc) (VectorXd, std::vector<double>&), 
	        problems choice, int* size);
