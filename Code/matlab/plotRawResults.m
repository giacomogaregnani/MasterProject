%% Load numerical solution Monte Carlo h variable N = 100
clear; clc; close all

for i = 0 : 3
    resultsfile{i+1} = ['../data/testHIRES_stableGAUSS_16_11_2016_04_38_', ...
        num2str(i), '.txt'];
end
nExperience = length(resultsfile);

%% Plot results MCMC

% numerical distribution h = 0.1
meanMCMC = zeros(nExperience, 5);
varMCMC = zeros(nExperience, 5);

for k = 1 : nExperience
    figure
    hold on
    results = dlmread(resultsfile{k});
    x = results(5001:end, :);
    meanMCMC(k, :) = mean(x);
    varMCMC(k, :) = var(x);
    plotmatrix(x)
end
set(gca, 'YTickLabel', '')
