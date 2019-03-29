function l = ProbLikelihood(w, xObs, uObs, obsNoise, f, x, nMC)

h = x(2) - x(1);

lVec = zeros(nMC, 1);

for i = 1 : nMC
    xProb = x;
    xProb(2:end-1) = xProb(2:end-1) + h^2 * (-0.5 + randn(size(xProb(2:end-1))));
    
    A = assembleMatrixProb(xProb);
    F = assembleRHSProb(f, xProb);
    
    uVec = [0; 1 / exp(w) * (A \ F); 0];
    u = interp1(x, uVec, xObs);
    
    diff = uObs - u;
    lVec(i) = exp(-0.5 / (obsNoise^2) * (diff' * diff));
end

l = log(mean(lVec));

% [maxL, maxI] = min(lVec);
% lVec(maxI) = [];
% 
% l = maxL + log(1 + sum(exp(lVec - maxL))) - log(nMC);

