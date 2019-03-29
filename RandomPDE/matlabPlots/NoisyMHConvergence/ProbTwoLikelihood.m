function l = ProbTwoLikelihood(w, xObs, uObs, obsNoise, A, F, x, nMC)

h = x(2) - x(1);

lVec = zeros(nMC, 1);
uVecDet = [0; 1 / exp(w) * (A \ F); 0];

for i = 1 : nMC    
    uVec = uVecDet + h * randn(size(uVecDet));
    u = interp1(x, uVec, xObs);
    
    diff = uObs - u;
    lVec(i) = exp(-0.5 / (obsNoise^2) * (diff' * diff));
end

[maxL, maxI] = min(lVec);
lVec(maxI) = [];

l = maxL + log(1 + sum(exp(lVec - maxL))) - log(nMC);

