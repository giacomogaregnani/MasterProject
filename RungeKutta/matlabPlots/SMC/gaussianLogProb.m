function g = gaussianLogProb(x, sigma, m)

g = -0.5 / sigma^2 * (x - m)' * (x - m);

end