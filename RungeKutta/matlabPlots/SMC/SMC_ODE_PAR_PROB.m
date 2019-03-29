function [sample, theta, weights] = SMC_ODE_PAR_PROB(nParticles, x0, fParam, RK, RK2, t0, tObs, yObs, h, obsProb, ...
                                                     a, nParam)
%% SMC for state estimation of an ODE model
% yObs matrix with an observation per column (each column is an R^d vector)
% theta is a matrix with a value per row

nObs = length(tObs);
tObs = [t0, tObs];

% ========= Initializations =========
weights = 1 / nParticles * ones(nParticles, 1);
sample = repmat(x0', nParticles, 1);
theta = randn(nParticles, nParam);
meanTheta = mean(theta);
varTheta = cov(theta);
s = sqrt(1 - a^2);

sampleMid = sample;
sampleNew = sample;
thetaNew = theta;
g = zeros(nParticles, 1);

k = 1;

while k <= nObs
    
    display(['time = ' num2str(tObs(k))])
    
    % ========= Propagation =========
    
    % Shrink the parameters
    theta = a * theta + (1 - a) * repmat(meanTheta, nParticles, 1);
    
    % Compute state predictor using RK
    nSteps = 1 / h * (tObs(k+1) - tObs(k));
    parfor i = 1 : nParticles
        f = @(t, y) fParam(t, y, theta(i, :));
        sampleMid(i, :) = RK(f, nSteps, sample(i, :)', h)';
    end
    
    % ========= Survival of the fittest =========
    
    % Compute fitness weights
    parfor i = 1 : nParticles
        g(i) = weights(i) * obsProb(sampleMid(i, :)', yObs(:, k));
    end
    
    S = sum(g);
    
    if S == 0 || isnan(S)
        
        display('resampling')
        idx = datasample(1:nParticles, nParticles, 'Weights', weights);
        sample = sample(idx, :);
        theta = theta(idx, :);
                
    else
        
        g = g / S;
             
        idx = datasample(1:nParticles, nParticles, 'Weights', g);
        
        % Reshuffle
        sample = sample(idx, :);
        sampleMid = sampleMid(idx, :);
        theta = theta(idx, :);
        
        % ========== Innovation ===========
        
        L = chol(varTheta);
        for i = 1 : nParticles
            thetaNew(i, :) = (theta(i, :)' + s * L' * randn(nParam, 1))';
        end
        
        parfor i = 1 : nParticles
            f = @(t, y) fParam(t, y, thetaNew(i, :));
            sampleNew(i, :) = RK2(f, nSteps, sample(i, :)', h)';
        end
        
        % ========== Weight updating ===========
        
        for i = 1 : nParticles
            weights(i) = obsProb(sampleNew(i, :)', yObs(:, k)) / obsProb(sampleMid(i, :)', yObs(:, k));
        end
        weights = weights / sum(weights);
        
        % ========== Update theta mean and var ===========
        
        meanTheta = sum(bsxfun(@times, thetaNew, weights));
        varTheta = zeros(nParam, nParam);
        for i = 1 : nParticles
            diff = thetaNew(i, :) - meanTheta;
            varTheta = varTheta + weights(i) * diff' * diff;
        end
        
        sample = sampleNew;
        theta = thetaNew;
        k = k + 1;
    end
end

end

