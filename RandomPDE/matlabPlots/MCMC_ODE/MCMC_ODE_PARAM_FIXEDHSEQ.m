clc; clear; close all
%%


filename = {'ODE_IC_PROB2_TEST_h_1';
    'ODE_IC_PROB_TEST_h_01';
    'ODE_IC_PROB_TEST_h_001'};
col = colormap(parula(8));
close
col = col(1:2:end, :);

% nPoints = 100;
% xEval = linspace(-2.4, -1.0, nPoints);
% yEval = linspace(-1.9, -1.5, nPoints);
% [XX, YY] = meshgrid(xEval, yEval);
xEval = linspace(-1.1, -0.8, 1000);
N = 100;

h = {'1', '01'};
MC = {'1', '5', '10', '100'};

for i = 1 : 2
    figure
    hold on
    for j = 1 : 4
        
        %     densTot = 0;
        %
        %     resultsTot = [];
        %     averages = [];
        %
        %     for i = 1 : N
        %
        %         display(i)
        %
        %         results = dlmread([filename{j}, '_exp_', num2str(i-1), '.txt']);
        %
        % %         resultsTot = [resultsTot; results];
        % %         averages = [averages; mean(results)];
        %
        %         %     [n, d] = size(results);
        %         %     bandwidth = 50 * std(results) * ((4 / (n * (d + 2)))^(1 / (d + 4)));
        %         %     m = mean(results);
        %
        %         dens = ksdensity(results, xEval);
        %         densTot = dens + densTot;
        %         %
        %         %     plot(xEval, dens, 'color', 0.7 * ones(3, 1));
        %     end
        %
        %     p(j) = plot(xEval, densTot / N, 'color', col(j, :), 'linewidth', 1);
        %
        resultsAndrea = dlmread(['METODOANDREA', '_h_', h{i}, '_MC_', MC{j}, '.txt']);
        dens = ksdensity(resultsAndrea, xEval);
        
        p(j) = plot(xEval, dens, 'color', col(j, :), 'linewidth', 1)
        
        resultsGiacomo = dlmread(['METODOGIACOMO', '_h_', h{i}, '_MC_', MC{j}, '.txt']);
        dens = ksdensity(resultsGiacomo, xEval);
        
        plot(xEval, dens, '--', 'color', col(j, :), 'linewidth', 1)
        
    end
    
    resultsDet = dlmread(['METODODET', '_h_', h{i}, '_MC_', MC{3}, '.txt']);
    dens = ksdensity(resultsDet, xEval);
    
    plot(xEval, dens, ':', 'color', col(j, :), 'linewidth', 1)
    
    legend([p(1), p(2), p(3), p(4)], {'$M = 1$', '$M = 5$', '$M = 10$', '$M = 100$'}, 'interpreter', 'laTeX', 'location', 'NE');
    
end

