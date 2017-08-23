function [thetaAll, accRatio, S] = InverseProblemPDE1D_MCMC_Andrea(t0, posterior, nMCMC, RAM, sigma, noisy)
%% Sample from parameter distribution using MCMC algorithm (Pseudo-marginal for probabilistic method)
% [thetaAll, accRatio, S] = InverseProblemPDE1D_MCMC(mesh, prior, f, nMCMC, prob, xObs, uObs, obsNoise)
% [t0] initial guess
% [prior] prior function
% [likelihood] likelihood function
% [nMCMC] the number of MCMC samples
% [prob] boolean (probabilistic or deterministic method for the PDE)
% [mesh] FEM mesh
% [RAM] boolean (adaptive MCMC)
% [sigma] Gaussian proposal std dev (or initial std dev for RAM)

% Evaluate prior and compute likelihood for initial guess
oldPosterior = posterior(t0);

% Initialize chain
nParam = length(t0);
thetaAll = zeros(nParam, nMCMC + 1);
thetaAll(:, 1) = t0;

% Initialize RAM
if size(sigma, 1) == 1;
    S = RAMInit(sigma, t0);
else
    S = sigma;
end
accRatio = 0;

if nargin == 5
   noisy = false; 
end

% MCMC loop

for i = 1 : nMCMC
    
    if mod(i, 100) == 0
        display(['iteration ' num2str(i), ', current acc. ratio ' num2str(accRatio / i,'%.4f') ...
            ', sigma = ', num2str(max(diag(S)),'%.4f'), ', posterior = ', num2str(oldPosterior, '%8.2E')]);
    end
    
    % Generate new guess (Gaussian proposal)
    w = randn(nParam, 1);
    theta = thetaAll(:, i) + S * w;
    
    % Evaluate prior and likelihood on new guess
    newPosterior = posterior(theta);
    
    if noisy
       oldPosterior = posterior(thetaAll(:, i)); 
    end
    
    % Compute alpha
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
    if RAM
        S = RAMUpdate(w, exp(alpha), i, 0.234, S);
    end
end

accRatio = accRatio / nMCMC;

end
