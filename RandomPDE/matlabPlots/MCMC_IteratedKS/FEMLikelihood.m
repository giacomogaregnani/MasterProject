function l = FEMLikelihood(w, xObs, uObs, obsNoise, A, F, x)

uVec = [0; 1 / exp(w) * (A \ F); 0];
u = interp1(x, uVec, xObs);

diff = uObs - u;
l = -0.5 / (obsNoise^2) * (diff' * diff);
