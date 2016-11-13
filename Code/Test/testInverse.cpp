#include <Solver.hpp>

int main()
{
    int n = 10;
    MatrixXd A = MatrixXd::Zero(10, 10);

    for (int i = 0; i < 10; i++) {
        for (int j = i; j < 10; j++) {
            A(j, i) = 1 + i + j;
        }
    }

    MatrixXd B = A.inverse();
    MatrixXd C = triInv(A, n);

    std::cout << A
              << std::endl
              << std::endl
              << B
              << std::endl
              << std::endl
              << C
              << std::endl;

    return 0;
}
