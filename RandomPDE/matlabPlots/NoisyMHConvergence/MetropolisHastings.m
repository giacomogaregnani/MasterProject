function [samples, accRatio] = MetropolisHastings(t0, posterior, proposal, nMCMC, noisy)
% Sample from parameter distribution using MCMC algorithm (Pseudo-marginal for probabilistic method)
% 
% [samples, accRatio] = MetropolisHastings(t0, posterior, proposal, nMCMC, noisy)
%
% INPUTS
% [t0]        initial guess
% [posterior] function handle to the posterior
% [proposal]  function handle to the proposal distribution
% [nMCMC]     the number of MCMC samples
% [noisy]     for the pseudo-marginal algorithm, perform a noisy update
%
% OUTPUTS
% [samples]   the Markov chain
% [accRatio]  final acceptance ratio

% Evaluate prior and compute likelihood for initial guess
oldPosterior = posterior(t0);

% Initialize chain
nParam = length(t0);
samples = zeros(nParam, nMCMC + 1);
samples(:, 1) = t0;

accRatio = 0;

if nargin == 4
    noisy = false;
end

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
    
    % For a noisy algorithm (MCWM), re-evaluate on old guess
    if noisy
       oldPosterior = posterior(samples(:, i)); 
    end
        
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
