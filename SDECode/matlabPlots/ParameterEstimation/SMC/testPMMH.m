clc; clear; close all;
%% Generate data

aTruth = 1.0; aMin = -1; aMax = 3;
S = 0.5;
IC = 1;
T = 1;
N = 100;

epsilon = 0.05;

f = @(a, x) -a * x;
fEps = @(a, x) -(a * x - 1/epsilon * sin(x/epsilon));
g = @(x) sqrt(2*S) * x.^0;

[t, obs] = EM_Traj(IC, @(x) fEps(aTruth, x), g, T, N, 1);

% Perturb
sigmaNoise = 1e-2;
noisyObs = [IC, obs(2:end) + sigmaNoise * randn(1, N)];

%% Posterior model

M = 1000;
sigmaProposal = 4e-2;
likelihood = @(a) PFLikelihood(a, IC, f, g, M, T, noisyObs, sigmaNoise);
prior = @(a) -0.5 * a^2;
posterior = @(a) likelihood(a) + prior(a);
proposal = @(a) gaussianProposal(a, sigmaProposal);

%% Compute posterior

isMCMC = 1;

if isMCMC
    nMCMC = 1000;
    initGuess = 0;
    sample = metropolisHastings(posterior, proposal, initGuess, nMCMC);
   
    %% Plot results
    [dens, xDens] = ksdensity(sample);
    
    figure
    plot(xDens, dens)
    hold on
    plot([aTruth, aTruth], [0, max(dens)])
else
    nSample = 200;
    sample = linspace(aMin, aMax, nSample);
    dens = zeros(size(sample));
    for i = 1 : nSample
        dens(i) = exp(posterior(sample(i)) / 100);
    end
    plot(sample, dens)
    hold on
    plot([aTruth, aTruth], [0, max(dens)])
end