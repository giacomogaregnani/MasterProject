%% Load numerical solution Monte Carlo h variable N = 100
clear; clc; close all

for i = 0 : 1
    resultsfile{i+1} = ['../data/TestBRUSS_SQRTGAUSS_11_11_2016_03_59', ...
        num2str(i), '.txt'];
end
nExperience = length(resultsfile);

%% Plot results MonteCarlo h variable N = 100

% numerical distribution h = 0.1
meanMCMC = zeros(nExperience, 1);
varMCMC = zeros(nExperience, 1);

for i = 1 : 1
    figure
    hold on
    for k = 1 : nExperience
        results = dlmread(resultsfile{k});
        x = results(5001:end, i);
        meanMCMC(k, i) = mean(x);
        varMCMC(k, i) = var(x);
        [f, xi] = ksdensity(x);
        plot(xi, f)
    end
    legend('h = 0.1', 'h = 0.05', 'h = 0.025', 'h = 0.0125', 'h = 0.0065', 'h = 0.0031')
end
set(gca, 'YTickLabel', '')
title('Gauss')
