function U = KL(uMean, V, sqrtD, theta)

U = uMean + sum(bsxfun(@times, theta', bsxfun(@times, (diag(sqrtD))', V)), 2);
