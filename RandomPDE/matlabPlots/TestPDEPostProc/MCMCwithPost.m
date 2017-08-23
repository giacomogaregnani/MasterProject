function [thetaAll, accRatio] = MCMCwithPost(posterior, theta, nMCMC, alphaStar)

% Evaluate
oldPosterior = posterior(theta);
nParam = length(theta);

% Initialize chain
thetaAll = zeros(nParam, nMCMC + 1);
thetaAll(:, 1) = theta;

% Initialize RAM
sigma = 0.1;
S = RAMInit(sigma, theta);
accRatio = 0;

% MCMC loop
for i = 1 : nMCMC
    
    if mod(i, 100) == 0
        display(['iteration ' num2str(i), ', current acc. ratio ' num2str(accRatio / i)...
            ', sigma = ', num2str(S(1,1)), ', (log)posterior = ', num2str(oldPosterior)]);
    end
    
    % Generate new guess (Gaussian proposal)
    w = randn(nParam, 1);
    theta = thetaAll(:, i) + S * w;
    
    % Compute alpha
    newPosterior = posterior(theta);
    alpha = newPosterior - oldPosterior;
    alpha = min(0, alpha);
    
    % Update the chain
    if alpha > log(rand(1))
        thetaAll(:, i+1) = theta;
        oldPosterior = newPosterior;
        accRatio = accRatio + 1;
    else
        thetaAll(:, i+1) = thetaAll(:, i);
    end
    
    % Update RAM matrix
    S = RAMUpdate(w, exp(alpha), i, alphaStar, S);
    
end
accRatio = accRatio / nMCMC;
