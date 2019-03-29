function L = computeLikelihoodEE(theta, tObs, yObs, sigmaObs, h)

L = 1;

for i = 1 : length(tObs)
    y = 1 + h * theta;
    L = L * exp(-1 / (2 * sigmaObs^2) * (yObs(i) - y^(tObs / h))^2);    
end
