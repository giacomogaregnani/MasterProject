function [sample, weights] = SMC_ODE(nPar, x0, RK, RK2, t0, tObs, yObs, h, obsProb, sample, weights)
%% SMC for state estimation of an ODE model
% yObs matrix with an observation per column (each column is an R^d vector)

nObs = length(tObs);
tObs = [t0, tObs];

% Initialization
if nargin == 9
    weights = 1 / nPar * ones(1, nPar);
    sample = x0 * ones(length(x0), nPar);
end
sampleNew = sample;
control = sample;
g = zeros(1, nPar);

for k = 1 : nObs
    
    % ========= Propagation =========
    
    nSteps = 1 / h * (tObs(k+1) - tObs(k));
    parfor i = 1 : nPar
        sample(:, i) = RK(nSteps, sample(:, i), h);
        [~, y] = RK2([tObs(k+1), tObs(k)], sample(:, i));
        control(:, i) = y(end);
    end
    
    % ========= Survival of the fittest =========
    
    % Compute fitness weights
    for i = 1 : nPar
        g(i) = weights(i) * exp(obsProb(sample(:, i), yObs(k)));
    end
    g = g / sum(g);
    
    % Draw indices with replacement
    idx = datasample(1:nPar, nPar, 'Weights', g);
    
    % Reshuffle
    sample = sample(:, idx);
    
    % ========== Innovation ===========
    
    for i = 1 : nPar
        sigma = norm(sample(:, i) - control(:, i));
        v = sigma * randn(size(sample(:, i)));
        sampleNew(:, i) = sample(:, i) + v;
    end
    
    % ========== Weight updating ===========
    
    for i = 1 : nPar
       weights(i) = obsProb(sampleNew(:, i), yObs(k)) / obsProb(sample(:, i), yObs(k)); 
    end
    weights = weights / sum(weights);
    
    sample = sampleNew;
    
end

end

