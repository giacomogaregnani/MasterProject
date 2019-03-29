function prior = GaussianPrior(w, priorMean, priorCovariance)

diff = priorMean - w;
prior = -0.5 * diff' * (priorCovariance \ diff);