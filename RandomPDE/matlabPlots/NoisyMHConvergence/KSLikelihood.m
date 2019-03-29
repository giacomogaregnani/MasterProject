function l = KSLikelihood(w, xObs, uObs, obsNoise, A, F, x, m, sigma)

uVec = [0; 1 / exp(w) * (A \ F); 0];
u = interp1(x, uVec, xObs);

diff = uObs - u - m;
l = -0.5 / (obsNoise^2 + sigma^2) * (diff' * diff);
