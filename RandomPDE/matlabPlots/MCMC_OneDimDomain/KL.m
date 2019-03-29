function U = KL(uMean, V, sqrtD, theta)

U = uMean + sum(bsxfun(@times, theta', bsxfun(@times, sqrtD', V)), 2);
