%% Parse file names
clear; clc; close all

resultsfile = cell(2, 50);
M = 50;
nExp = 100;
nM = 1;

for i = 1 : nM
    for j = 1 : nExp
        resultsfile{i, j} = ['../data/testFITZNAG_MCvar_21_11_2016_02_47__MC_', num2str(M(i)), '_h_10_Exp_', num2str(j - 1), '.txt'];
    end
end

%% Compute means and variance of the mean

meanLik = zeros(nM, nExp);
varMC = zeros(nM, 1);
varLik = zeros(nM, 1);

for i = 1 : nM
    for j = 1 : nExp
        results = dlmread(resultsfile{i, j});
        x = results(:, 1 : 3);
        normX = sqrt(sum(x.^2, 2));
        likelihood = results(:, 4);
        okSamples = find(likelihood > -1e4);
        x = x(okSamples, : );
        likelihood = likelihood(okSamples);
%         if (j < 10)
%             figure
%             plot(x(:, 1), likelihood, 'o')
%         end
        meanMC = mean(x);
        meanLik(i, j) = mean(likelihood);
    end
    varLik(i) = var(meanLik(i, :));
end

[f, xi] = ksdensity(meanLik - mean(meanLik));
plot(xi, f)
