function L = computeLikelihood(theta, tObs, yObs, sigmaObs)

L = 1;

for i = 1 : length(tObs)
    L = L * exp(-1 / (2 * sigmaObs^2) * (yObs(i) - exp(theta * tObs(i)))^2);    
end
