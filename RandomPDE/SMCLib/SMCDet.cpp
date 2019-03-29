#include "SMCDet.hpp"

void SMCDet::updateMeanAndCovariance(std::vector<VectorXd>& vec,
                                  std::vector<double>& weights,
                                  VectorXd& mean, MatrixXd& cov)
{
    mean.setZero();
    for (unsigned long i = 0; i < nParticles; i++) {
        mean += weights[i] * vec[i];
    }

    cov.setZero();
    VectorXd tmp;
    for (unsigned long i = 0; i < nParticles; i++) {
        tmp = vec[i] - mean;
        cov += weights[i] * tmp * tmp.transpose();
    }
}

double GaussianPDFdet(VectorXd& x, VectorXd& m, double sigma)
{
    return std::exp(-0.5 / (sigma * sigma) * (x - m).dot(x - m));
}


SMCDet::SMCDet(unsigned long N, double a,
               VectorXd &IC, odeDef ODE, Butcher tableau,
               Butcher tableau2, double h,
               std::default_random_engine* generator,
               VectorXd &priorMean, double priorStdDev):
        nParticles(N),
        a(a),
        priorMean(priorMean),
        priorStdDev(priorStdDev),
        IC(IC),
        generator(generator),
        h(h)
{
    detMethod = RungeKutta(ODE, tableau);
    detMethod2 = RungeKutta(ODE, tableau2);
    sizeParam = priorMean.size();
}

void SMCDet::compute(std::vector<VectorXd>& samples,
                  std::vector<VectorXd>& thetas,
                  std::vector<double>& weights,
                  std::vector<double>& tObs,
                  std::vector<VectorXd>& yObs,
                  double obsNoise)
{
    // Initialization
    samples = std::vector<VectorXd>(nParticles, IC);
    weights = std::vector<double>(nParticles, 1.0 / nParticles);
    thetas = std::vector<VectorXd>(nParticles, priorMean);

    std::vector<VectorXd> midSamples(nParticles, IC);
    std::vector<VectorXd> newSamples(nParticles, IC);
    std::vector<VectorXd> newThetas(nParticles, priorMean);
    std::vector<double> g(nParticles, 0.0);
    std::vector<unsigned long> replacementIdx(nParticles);
    std::map<int, int> replacementMap;

    for (unsigned long i = 0; i < nParticles; i++) {
        for (int j = 0; j < sizeParam; j++) {
            thetas[i](j) += priorStdDev * gaussian(*generator);
        }
    }
    VectorXd meanTheta(sizeParam);
    MatrixXd covTheta(sizeParam, sizeParam);
    updateMeanAndCovariance(thetas, weights, meanTheta, covTheta);
    double s = std::sqrt(1 - a * a);

    VectorXd control(IC.size()), randomContr(IC.size()), sigma(IC.size());

    tObs.insert(tObs.begin(), 0.0);
    unsigned long nSteps;

    for (unsigned int k = 0; k < yObs.size(); k++) {

        std::cout << "time = " << tObs[k] << std::endl;

        // STEP 1 : PROPAGATION
        for (unsigned long i = 0; i < nParticles; i++) {
            thetas[i] = a * thetas[i] + (1 - a) * meanTheta;
        }

        nSteps = static_cast<unsigned long>(std::round((tObs[k+1] - tObs[k]) / h));
        VectorXd y(samples[0].size());
        for (unsigned long j = 0; j < nParticles; j++) {
            y = samples[j];
            for (unsigned long i = 0; i < nSteps; i++) {
                y = detMethod.oneStep(h, y, thetas[j]);
            }
            midSamples[j] = y;
        }

        // STEP 2 : SURVIVAL OF THE FITTEST
        double sum = 0.0;
        for (unsigned long i = 0; i < nParticles; i++) {
            g[i] = weights[i] * GaussianPDFdet(midSamples[i], yObs[k], obsNoise);
            sum += g[i];
        }
        for (unsigned long i = 0; i < nParticles; i++) {
            g[i] /= sum;
        }

        // shuffle
        std::discrete_distribution<unsigned long> replacementGenerator(g.begin(), g.end());
        for(int i = 0; i < nParticles; i++) {
            unsigned long idx = replacementGenerator(*generator);
            samples[i] = samples[idx];
            thetas[i] = thetas[idx];
            midSamples[i] = midSamples[idx];
        }

        // STEP 3 : INNOVATION
        LLT<MatrixXd> covCholesky(covTheta);
        MatrixXd L = covCholesky.matrixL();
        for (unsigned long i = 0; i < nParticles; i++) {
            VectorXd randomVec(sizeParam);
            for (int j = 0; j < sizeParam; j++) {
                randomVec(j) = priorStdDev * gaussian(*generator);
            }
            newThetas[i] = thetas[i] + s * (L * randomVec);
        }

        for (unsigned long j = 0; j < nParticles; j++) {

            control = samples[j];
            for (unsigned long i = 0; i < nSteps; i++) {
                control = detMethod2.oneStep(h, control, thetas[j]);
            }

            for (int indSize = 0; indSize < control.size(); indSize++) {
                sigma(indSize) = (control(indSize) - midSamples[j](indSize));
                randomContr(indSize) = std::abs(sigma(indSize)) * gaussian(*generator);
            }
            // std::cout << randomContr.transpose() << std::endl;

            newSamples[j] = midSamples[j] + randomContr;
        }

        // STEP 4 : WEIGHT UPDATING
        sum = 0.0;
        for (unsigned long j = 0; j < nParticles; j++) {
            weights[j] = GaussianPDFdet(newSamples[j], yObs[k], obsNoise) / GaussianPDFdet(midSamples[j], yObs[k], obsNoise);
            sum += weights[j];
        }
        for (unsigned long j = 0; j < nParticles; j++) {
            weights[j] /= sum;
        }

        // STEP 5 : UPDATE MEAN AND VAR
        updateMeanAndCovariance(newThetas, weights, meanTheta, covTheta);

        samples = newSamples;
        thetas = newThetas;
    }
    tObs.erase(tObs.begin());
}
