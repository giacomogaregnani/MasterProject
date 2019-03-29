function [m, sigma] = estimateErrorStats(nTrials, xObs, A, F, x, priorMean, priorVariance, uEx)

err = zeros(nTrials, 1);
priorStdDev = sqrt(priorVariance);

for i = 1 : nTrials
    
    w = priorMean + priorStdDev * randn(1);
    
    uVec = [0; 1 / exp(w) * (A \ F); 0];
    uFEM = interp1(x, uVec, xObs);
    u = 1 / exp(w) * uEx(xObs);
    
    err(i) = u - uFEM;    
    
end

m = mean(err);
sigma = var(err);