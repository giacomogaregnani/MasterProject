#include <KarhunenLoeve.hpp>
#include <fstream>

int main(void)
{
    int N = 100;

    VectorXd theta(2);
    theta << 0.0, 1.0;

    VectorXd KLMean = VectorXd::Ones(N);
    KarhunenLoeve KL(KLMean, INVLAPDIR, N);
    VectorXd field = KL.KL(theta);

    std::ofstream output(DATA_PATH + std::string("TestKL.txt"), std::ofstream::out | std::ofstream::trunc);
    output << field.transpose() << std::endl;
    output.close();

    return 0;
}
