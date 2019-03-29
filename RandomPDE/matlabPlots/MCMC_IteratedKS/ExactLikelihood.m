function l = ExactLikelihood(w, xObs, uObs, obsNoise)

u = 1 / (exp(w) * (2 * pi)^2) * sin(2 * pi * xObs);

l = -0.5 / (obsNoise^2) * (uObs - u)' * (uObs - u);
