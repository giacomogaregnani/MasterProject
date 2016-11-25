%% Load numerical solution Monte Carlo h variable N = 100
clear; clc; close all

resultsfile = cell(3, 3);
h = [0.2, 0.1, 0.05, 0.025];
M = [1, 10, 100, 1000];
hText = {'0_2', '0_1', '0_05', '0_025'};
[nH, nM] = size(resultsfile);

for i = 1 : nH
    for j = 1 : nM
        resultsfile{i,j} = ['../data/testBRUSS_MC_19_11_2016_03_02__MC_', num2str(M(j)), '_h_', hText{i}, '.txt'];
    end
end

%% Plot distributions

meanMC = zeros(nH, nM);
errors = ones(nH, nM);

for i = 1 : nH
    figure
    hold on
    leg = cell(nM, 1);
    for j = 1 : nM
        results = dlmread(resultsfile{i, j});
        x = results(5001:end);
        meanMC(i, j) = mean(x);
        [f, xi] = ksdensity(x);
        plot(xi, f)
        leg{j} = ['M = ', num2str(M(j))];
        errors(i, j) = abs(errors(i, j) - meanMC(i, j));
    end
    legend(leg);
end
set(gca, 'YTickLabel', '')









