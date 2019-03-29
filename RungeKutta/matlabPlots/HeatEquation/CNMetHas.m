function [thetaAll, accRatio] = CNMetHas(t0, posterior, nMCMC, proposal, noisy)
% Giacomo Garegnani (EPFL SB MATH ANMC) 2017
% Reference: MCMC Methods for Functions: Modifying old algorithms to make
% them faster, Cotter, Roberts, Stuart and White, 2013
%
% Sample from parameter distribution using Crank Nicolson Metropolis algorithm 
% 
% [thetaAll, accRatio] = CNMetHas(t0, posterior, nMCMC, proposal, noisy)
%
% INPUTS
% [t0]        initial guess
% [posterior] function handle to the posterior
% [nMCMC]     the number of MCMC samples
% [proposal]  anonymous function for proposal distribution
% [noisy]     booelan: perform a noisy update (pseudo-marginal)
%
% OUTPUTS
% [thetaAll]  the Markov chain
% [accRatio]  final acceptance ratio

% Evaluate prior and compute likelihood for initial guess
oldPosterior = posterior(t0);

% Initialize chain
nParam = length(t0);
thetaAll = zeros(nParam, nMCMC + 1);
thetaAll(:, 1) = t0;

% Initialize RAM
accRatio = 0;

if nargin == 4
   noisy = false; 
end

% MCMC loop

for i = 1 : nMCMC
    
    if mod(i, 100) == 0
        display(['iteration ' num2str(i), ', current acc. ratio ' num2str(accRatio / i,'%.4f') ...
            ', posterior = ', num2str(oldPosterior, '%8.2E')]);
    end
    
    % Generate new guess
    theta = proposal(thetaAll(:, i));
    
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
end

accRatio = accRatio / nMCMC;

end
