function g = gaussianProb(x, sigma, m)

g = exp(-0.5 / sigma^2 * (x - m)' * (x - m));

end