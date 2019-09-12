#include "NoiselessPosterior.hpp"

NLPosterior::NLPosterior(std::vector<double> &x, double h,
                         double (*gradV0) (double), VectorXd& paramSde):
        obs(x),
        h(h),
        gradV0(gradV0),
        paramSde(paramSde)
{
    N = obs.size() - 1;
}

double NLPosterior::computePosterior(VectorXd &theta)
{
    double A = theta(0);
    double ASqd = A * A;
    double prior = -0.5 * ASqd;

    paramSde(1) = A;

    double likelihood = 0.0;
    for (unsigned long i = 0; i < N; i++) {
        double gradLoc = gradV0(obs[i]);
        likelihood -= A * gradLoc * (obs[i+1] - obs[i])
                      + 0.5 * h * ASqd * gradLoc * gradLoc;
    }

    double posterior = prior + likelihood;
    return posterior;
}