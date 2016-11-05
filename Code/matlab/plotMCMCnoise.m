% MCMC results
%    results = dlmread('../data/resultsNoise.txt');


clc
clear
% close all

% resultsfile = {'resultsMCMC_Noise_h01.mat', ...
%     'resultsMCMC_Noise_h005.mat', ...
%     'resultsMCMC_Noise_h002.mat', ...
%     'resultsMCMC_Noise_h001.mat', ...
%     'resultsMCMC_Noise_h0005.mat'
%     };

resultsfile = {'resultsMCMC_Noise_h01.mat', ...
    'resultsMCMC_Noise_h005.mat', ...
    'resultsMCMC_Noise_h002.mat', ...
    'resultsMCMC_Noise_h001.mat', ...
    'resultsMCMC_Noise_h0005.mat'
    };

col = {'b', 'r', 'k', 'm', 'c'};
meanMCMC = zeros(length(resultsfile), 3);
varMCMC = zeros(length(resultsfile), 3);
figure
for k = 1 : 5
    load(resultsfile{k});
    x = results(10001:end, :);
    meanMCMC(k, :) = mean(x);
    varMCMC(k, :) = var(x);
    [n, d] = size(x);
    for j = 1 : d
        for i = 1 : j
            subplot(d, d, d * (j - 1) + i);
            hold on
            if i == j
                %                 h = histogram(x(:, j), 'Normalization', 'pdf', 'DisplayStyle', 'stairs');
                [f, xi] = ksdensity(x(:, j));
                plot(xi, f, col{k})
                set(gca, 'YTickLabel', '')
            else
                [first,second] = hist3([x(:, i), x(:, j)]);
                contour(second{1},second{2},first,col{k})
                if j == 2
                    ylim([-0.25 0.4])
                end
                if j == 3
                    ylim([2.5 3.1])
                end
            end
            if i == 1
                xlim([0.05 0.3])
            end
            if i == 2
                xlim([-0.25 0.4])
            end
            if i == 3
                xlim([2.5 3.1])
            end
        end
    end
end
hL = legend('h = 0.1', 'h = 0.05', 'h = 0.02', 'h = 0.01', 'h = 0.005');
newPosition = [0.75 0.75 0.2 0.2];
newUnits = 'normalized';
set(hL, 'Position', newPosition, 'Units', newUnits);

%% all plots for report
% MCMC results
%    results = dlmread('../data/resultsNoise.txt');


clc
clear
% close all

resultsfile = {'resultsMCMC_Noise_h01.mat', ...
    'resultsMCMC_Noise_h005.mat', ...
    'resultsMCMC_Noise_h002.mat', ...
    'resultsMCMC_Noise_h001.mat', ...
    'resultsMCMC_Noise_h0005.mat'
    };

col = {'b', 'r', 'k', 'm', 'c'};
meanMCMC = zeros(length(resultsfile), 3);
varMCMC = zeros(length(resultsfile), 3);

load(resultsfile{1})

[n, d] = size(results(10001:end, :));
for j = 1 : d
    for i = 1 : j
        figure
        for k = 1 : 5
            load(resultsfile{k});
            x = results(10001:end, :);
            if i == j
                [f, xi] = ksdensity(x(:, j));
                pEq = plot(xi, f, col{k});
                hold on
                set(gca, 'YTickLabel', '')
                if i == 1 && k == 5
                   legend('h = 0.1', 'h = 0.05', 'h = 0.02', 'h = 0.01', 'h = 0.005');
                end
            else
                [first,second] = hist3([x(:, i), x(:, j)]);
                pDiff = contour(second{1},second{2},first,col{k});
                hold on
                if j == 2
                    ylim([-0.25 0.4])
                end
                if j == 3
                    ylim([2.5 3.1])
                end
            end
            if i == 1
                xlim([0.05 0.3])
            end
            if i == 2
                xlim([-0.25 0.4])
            end
            if i == 3
                xlim([2.5 3.1])
            end
        end
    end
end


