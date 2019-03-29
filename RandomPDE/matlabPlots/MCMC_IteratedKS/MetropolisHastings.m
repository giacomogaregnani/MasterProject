function [samples, accRatio] = MetropolisHastings(t0, posterior, proposal, nMCMC)
% Sample from parameter distribution using MCMC algorithm (Pseudo-marginal for probabilistic method)
% 
% [thetaAll, accRatio] = InverseProblemPDE1D_MCMC(t0, posterior, proposal, nMCMC)
%
% INPUTS
% [t0]        initial guess
% [posterior] function handle to the posterior
% [proposal]  proposal distribution
% [nMCMC]     the number of MCMC samples
%
% OUTPUTS
% [thetaAll]  the Markov chain
% [accRatio]  final acceptance ratio

% Evaluate prior and compute likelihood for initial guess
oldPosterior = posterior(t0);

% Initialize chain
nParam = length(t0);
samples = zeros(nParam, nMCMC + 1);
samples(:, 1) = t0;

accRatio = 0;

% MCMC loop
for i = 1 : nMCMC
    
    if mod(i, 1000) == 0
        display(['iteration ' num2str(i), ', current acc. ratio ' num2str(accRatio / i,'%.4f') ...
                 ', posterior = ', num2str(oldPosterior, '%8.2E')]);
    end
    
    % Generate new guess (Gaussian proposal)
    proposed = proposal(samples(:, i));
    
    % Evaluate posterior on new guess
    newPosterior = posterior(proposed);
        
    % Compute alpha
    alpha = newPosterior - oldPosterior;
    alpha = min(0, alpha);
    
    % Update the chain
    
    if alpha > log(rand(1))
        samples(:, i+1) = proposed;
        oldPosterior = newPosterior;
        accRatio = accRatio + 1;
    else
        samples(:, i+1) = samples(:, i);
    end
   
end

accRatio = accRatio / nMCMC;

end
