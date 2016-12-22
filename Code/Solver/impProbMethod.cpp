#include "Solver.hpp"

impProbMethod::impProbMethod(odeDef ODE, double timestep, std::vector<double> paramVec,
                             double stoch, MatrixXd A, VectorXd b, int order)
{
    detSolver = std::make_shared<ImplicitRK>(ODE, paramVec, A, b, order);
    solution = ODE.initialCond;
    sigma = stoch;
    rootsigma = sqrt(sigma);
    size = ODE.size;
    h = timestep;
    int detOrder = detSolver->getOrder();
    hfunc = pow(h, detOrder + 0.5);
    IC = ODE.initialCond;
    parameters = paramVec;
}

VectorXd& impProbMethod::getSolution(void)
{
    return solution;
}

void impProbMethod::oneStep(std::default_random_engine& generator, double step)
{
    VectorXd noise(size);
    for (int i = 0; i < size; i++){
        noise(i) = normalDist(generator);
    }

    if (step == h) {
        VectorXd detSol = detSolver->oneStep(solution, h);
        solution = detSol + hfunc * rootsigma * noise;
    } else {
        double hfunctmp = pow(h, detSolver->getOrder() + 0.5);
        VectorXd detSol = detSolver->oneStep(solution, step);
        solution = detSol + hfunctmp * rootsigma * noise;
    }
}

VectorXd impProbMethod::oneStepGiven(VectorXd& noise, double step)
{
    if (step == h) {
        VectorXd detSol = detSolver->oneStep(solution, h);
        solution = detSol + hfunc * rootsigma * noise;
    } else {
        double hfunctmp = pow(h, detSolver->getOrder() + 0.5);
        VectorXd detSol = detSolver->oneStep(solution, step);
        solution = detSol + hfunctmp * rootsigma * noise;
    }
    return solution;
}

void impProbMethod::resetIC(void)
{
    solution = IC;
}
