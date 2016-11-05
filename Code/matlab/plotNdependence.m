% MCMC results
clc
clear
close all

%% Important dataset
% ../data/testMCone (0 - 5) MC on one data at T = 1 for h = 0.1
% ../data/testML (0 - 3) MLMC on one data at T = 1
%   for h = {0.1, 0.05, 0.02, 0.01}

%%

% resultsfile{1} = ['../data/test', ...
%         num2str(0), '.txt'];
for i = 0 : 4
    resultsfile{i+1} = ['../data/testRAM_GAUSS31_10_2016_04_18', ...
        num2str(i), '.txt'];
end
% resultsfile{1} = ['../data/testRAM_GAUSS_T1_h00131_10_2016_04_09', ...
%          num2str(0), '.txt'];
    
% small set of data
% for i = 0 : 3
%     resultsfile{i+1} = ['../data/testML', ...
%         num2str(i), '.txt'];
% end

% resultsfile{5} = ['../data/testMCone', ...
%         num2str(5), '.txt'];

% H = 0.005
% for i = 0 : 3
%     resultsfile{i+1} = ['../data/solNMC002H', ...
%         num2str(i), '.txt'];
% end

% H = 0.02
% for i = 0 : 4
%     resultsfile{i+1} = ['../data/solNMC', ...
%         num2str(i), '.txt'];
% end

% H = 0.1
% for i = 0 : 4
%     resultsfile{i+1} = ['../data/solNMC01H', ...
%         num2str(i), '.txt'];
% end

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
%                 if j == 2
%                     ylim([-0.25 0.4])
%                 end
%                 if j == 3
%                     ylim([2.5 3.1])
%                 end
            end
%             if i == 1
%                 xlim([0.05 0.3])
%             end
%             if i == 2
%                 xlim([-0.25 0.4])
%             end
%             if i == 3
%                 xlim([2.5 3.1])
%             end
        end
    end
end
hL = legend('N = 1', 'N = 10', 'N = 100', 'N = 1000', 'N = 10000', 'N = 100000');
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


