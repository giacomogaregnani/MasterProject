function l = PFLikelihood(a, IC, f, g, M, T, noisyObs, sigmaNoise)

[x, w, l] = PF(IC, @(x)f(a, x), g, M, T, noisyObs, sigmaNoise);

% Choose the most likely path
idx = randsample(1:M, 1, true, w);
xChosen = x(idx, :);

t = linspace(0, 1, length(xChosen));
plot(t, xChosen)
hold on
plot(t, noisyObs)
return