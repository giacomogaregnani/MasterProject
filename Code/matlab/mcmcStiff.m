% MCMC results
clc
clear
close all

for i = 0 : 3
   resultsfile{i+1} = ['../data/solNMC', ...
                       num2str(i), '.txt'];
end
       
figure
hold on
col = {'b', 'r', 'k', 'm', 'c'};       
meanMCMC = zeros(length(resultsfile), 3);
varMCMC = zeros(length(resultsfile), 3);


for k = 1 : length(resultsfile);
    results = dlmread(resultsfile{k});
    x = results(10001:end, :);
    meanMCMC(k, :) = mean(x);
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
                [first,second] = hist3([x(:, i), x(:, j)]);
                contour(second{1}, second{2}, first, col{k})
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
%%
% for i = 0 : 2
%    resultsfilestiff{i+1} = ['../data/resultsLorenzStab', ...
%                        num2str(i), '.txt'];
% end
%        
% figure
% col = {'b', 'r', 'k', 'm', 'c'};       
% meanMCMC = zeros(length(resultsfilestiff), 3);
% varMCMC = zeros(length(resultsfilestiff), 3);
% for k = 1 : length(resultsfilestiff);
%     results = dlmread(resultsfilestiff{k});
%     x = results(10001:end, :);
%     meanMCMC(k, :) = mean(x);
%     varMCMC(k, :) = var(x);
%     [n, d] = size(x);
%     for j = 1 : d
%         for i = 1 : j
%             subplot(d, d, d * (j - 1) + i);
%             hold on
%             if i == j
% %                 h = histogram(x(:, j), 'Normalization', 'pdf', 'DisplayStyle', 'stairs');
%                 [f, xi] = ksdensity(x(:, j));
%                 plot(xi, f, col{k})
%                 set(gca, 'YTickLabel', '')
%             else
%                 [first,second] = hist3([x(:, i), x(:, j)]);
%                 contour(second{1},second{2},first,col{k})
%             end
%         end
%     end
% end
% hL = legend('h = 0.1', 'h = 0.05', 'h = 0.02', 'h = 0.01', 'h = 0.005');
% newPosition = [0.75 0.75 0.2 0.2];
% newUnits = 'normalized';
% set(hL, 'Position', newPosition, 'Units', newUnits);
% 
% 
% 
% 
