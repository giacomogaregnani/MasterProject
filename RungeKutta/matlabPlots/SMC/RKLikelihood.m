function l = RKLikelihood(y, yObs, tObs, fParam, RK, h, obsProb, w)

l = 0;
tObs = [0, tObs];

for k = 1 : length(tObs) - 1
    
    nSteps = 1 / h * (tObs(k+1) - tObs(k));
    f = @(t, y) fParam(t, y, w);
    
    for i = 1 : nSteps
        y = RK(f, nSteps, y, h);
    end
    
    l = l + obsProb(y, yObs(k));
end
