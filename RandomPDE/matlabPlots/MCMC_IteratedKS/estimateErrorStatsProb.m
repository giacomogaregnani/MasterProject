function [m, sigma] = estimateErrorStatsProb(nTrials, xObs, x, f, priorMean, priorVariance, uEx, nMC)

err = zeros(nTrials, 1);
varVector = zeros(nTrials, 1);
priorStdDev = sqrt(priorVariance);
h = x(2) - x(1);

for i = 1 : nTrials
    
    w = priorMean + priorStdDev * randn(1);
    
    uFEM = zeros(nMC, length(xObs));
       
    for j = 1 : nMC    
        xProb = x;
        xProb(2:end-1) = xProb(2:end-1) + h^2 / 2 * (-0.5 + randn(size(xProb(2:end-1))));
        
        A = assembleMatrixProb(xProb);
        F = assembleRHSProb(f, xProb);
        
        uVec = [0; 1 / exp(w) * (A \ F); 0];
        
        uFEM(j, :) = interp1(x, uVec, xObs);        
    end
  
    u = 1 / exp(w) * uEx(xObs);
    err(i) = sqrt(var(uFEM));
  
end

m = mean(err);
sigma = var(err);