function [sample, theta, weights] = SMC_ODE_PAR(nParticles, x0, fParam, RK, RK2, t0, tObs, yObs, h, obsProb, ...
                                                a, nParam)
%% SMC for state estimation of an ODE model
% yObs matrix with an observation per column (each column is an R^d vector)
% theta is a matrix with a value per row

nObs = length(tObs);
tObs = [t0, tObs];

% Initialization
weights = 1 / nParticles * ones(nParticles, 1);
sample = x0 * ones(nParticles, length(x0));
theta = randn(nParticles, nParam);
meanTheta = mean(theta);
varTheta = cov(theta);
s = sqrt(1 - a^2);

sampleMid = sample;
sampleNew = sample;
control = sample;
thetaNew = theta;
g = zeros(nParticles, 1);

for k = 1 : nObs
    
    display(['time = ' num2str(tObs(k))])
    
    % ========= Propagation =========
    
    % Shrink the parameters
    theta = a * theta + (1 - a) * repmat(meanTheta, nParticles, 1);
    nSteps = 1 / h * (tObs(k+1) - tObs(k));
    
    parfor i = 1 : nParticles
        f = @(t, y) fParam(t, y, theta(i, :));
        sampleMid(i, :) = RK(f, nSteps, sample(i, :), h);
    end
    
    % ========= Survival of the fittest =========
    
    % Compute fitness weights
    parfor i = 1 : nParticles
        g(i) = weights(i) * obsProb(sampleMid(i, :), yObs(k));
    end
    g = g / sum(g);
    
    % Draw indices with replacement
    idx = datasample(1:nParticles, nParticles, 'Weights', g);
    
    % Reshuffle
    sampleMid = sampleMid(idx, :);
    theta = theta(idx, :);
    
    % ========== Innovation ===========
    
    L = chol(varTheta);
    for i = 1 : nParticles
        thetaNew(i, :) = (theta(i, :)' + s * L' * randn(nParam, 1))';
    end
    
    parfor i = 1 : nParticles
        f = @(t, y) fParam(t, y, thetaNew(i, :));
        control(i, :) = RK2(f, nSteps, sample(i, :), h);
        sigma = (sampleMid(i, :) - control(i, :)).^2;
        v = (diag(sigma) * randn(size(sample(i, :)))')';
        sampleNew(i, :) = sampleMid(i, :) + v;
    end
    % ========== Weight updating ===========
    
    for i = 1 : nParticles
        weights(i) = obsProb(sampleNew(i, :), yObs(k)) / obsProb(sampleMid(i, :), yObs(k));
    end
    weights = weights / sum(weights);
    
    % ========== Update theta mean and var ===========
    meanTheta = sum(bsxfun(@times, thetaNew, weights));
    varTheta = zeros(nParam, nParam);
    for i = 1 : nParticles
        varTheta = varTheta + weights(i) * thetaNew(i, :)' * thetaNew(i, :);
    end
    
    sample = sampleNew;
    theta = thetaNew;
end

end

