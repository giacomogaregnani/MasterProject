function l = likelihoodProbFEM(param, field, f, x, xObs, observations, sigma, M)

l = zeros(M, 1);
p = 1;

paramField = @(xx) field(xx, param);

parfor i = 1 : M
    
    X = x;

    for kk = 2 : length(x)-1
        hBar = min(x(kk)-x(kk-1), x(kk+1)-x(kk));
        X(kk) = x(kk) + hBar^p * (rand - 0.5);
    end
    
    F = assembleRHS(f, X);
    A = assembleMatrix(paramField, X);
    uVec = [0; A \ F; 0];
    
    u = interp1(X, uVec, xObs);
    
    diff = observations - u;
    l(i) = -1 / (2 * sigma^2) * (diff') * diff;
    
end

[maxLik, maxI] = max(l); 
l = l([1:maxI-1,maxI+1:end]);
others = sum(exp(l - maxLik));
l = maxLik + log(1 + others) - log(M);

end