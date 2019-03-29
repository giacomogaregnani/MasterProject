function likelihood = likelihood_CN(u0, A, tSpan, h, xObs, yObs, x, C)

y = impEulerLinear(A, u0, tSpan, h);

uPoints = interp1(x(2:end-1), y, xObs);

% Compute log-likelihood
likelihood = -0.5 * (uPoints - yObs) * (C \ (uPoints - yObs))';

