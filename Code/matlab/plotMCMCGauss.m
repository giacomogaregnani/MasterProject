% MCMC results
clc
clear
close all

%% Important dataset
% ../data/testMCone (0 - 5) MC on one data at T = 1 for h = 0.1
% ../data/testML (0 - 3) MLMC on one data at T = 1
%   for h = {0.1, 0.05, 0.02, 0.01}

%%
for i = 0 : 2
    resultsfile{i+1} = ['../data/testRAM_GAUSS31_10_2016_04_47', ...
        num2str(i), '.txt'];
end

figure
hold on
col = {'b', 'r', 'k', 'm', 'c', 'g'};
meanMCMC = zeros(length(resultsfile), 3);
varMCMC = zeros(length(resultsfile), 3);

for k = 1 : length(resultsfile);
    results = dlmread(resultsfile{k});
    x = results(10001:end, :);
    meanMCMC(k, :) = mean(x);
    C{k} = cov(x);
    varMCMC(k, :) = var(x);
    [n, d] = size(x);
    for j = 1 : d
        for i = 1 : j
            subplot(d, d, d * (j - 1) + i);
            hold on
            if i == j
                [f, xi] = ksdensity(x(:, j));
                plot(xi, f, col{k})
                set(gca, 'YTickLabel', '')
            else
%                 scatter(x(:, i), x(:, j), col{k})
                %                 binScatterPlot(x(:, i), x(:, j), col{k});
%                 [first,second] = hist3([x(:, i), x(:, j)]);
%                 contour(second{1}, second{2}, first, col{k})
            end
        end
    end
end
hL = legend('h = 0.05', 'h = 0.02', 'h = 0.01', 'h = 0.005');
newPosition = [0.75 0.75 0.2 0.2];
newUnits = 'normalized';
set(hL, 'Position', newPosition, 'Units', newUnits);


%% Hellinger distance

nEx = length(resultsfile);
detF = det(C{nEx});
H = zeros(1, length(resultsfile)-2);
for i = 1 : nEx - 1
    m = meanMCMC(i, :) - meanMCMC(end, :);
    MATRIX = 0.5 * (C{i} + C{end});
    H(i) = 1 - (det(C{i})^0.25 * detF^0.25 / (det(MATRIX)^0.5) * ...
        exp(-1/8 * m * (MATRIX \ m')));
    H(i) = sqrt(H(i));
end

figure
N = 10.^[0:length(resultsfile)-2];
loglog(N, H,'o-')
hold on
delta = 0.0;
loglog(N, log(N) ./ (N.^(0.5 - delta)),'o-')


