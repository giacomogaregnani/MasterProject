clc; clear; close all

truth = [0.5, 0];
col = colormap(parula(10));
col = col(1:2:end, :);
close
enhanced = 1;
fontsizeLAB = getLatexTextSize('normalsize', 'enhanced', enhanced);
fontsizeTICK = getLatexTextSize('footnotesize', 'enhanced', enhanced);
W = 7; H = 7;

createFigure(W, H, 'enhanced',enhanced);
hold on
i = 1;
for fileName = {'ODE_IC_DET_KEPLER_h_02.txt', 'ODE_IC_DET_KEPLER_h_01.txt', 'ODE_IC_DET_KEPLER_h_005.txt', 'ODE_IC_DET_KEPLER_h_0025.txt', 'ODE_IC_DET_KEPLER_h_small.txt'}

    results = dlmread(fileName{1});
   
    % Evaluate 2D density
    nPoints = 50;
    xEval = linspace(min(results(:, 3)), max(results(:, 3)), nPoints);
    yEval = linspace(min(results(:, 4)), max(results(:, 4)), nPoints);
    [XX, YY] = meshgrid(xEval, yEval);
    
    [n, d] = size(results);    
    bandwidth = std(results) * ((4 / (n * (d + 2)))^(1 / (d + 4)));
    m = mean(results);
        
    dens = mvksdensity(results, [m(1) * ones(size(XX(:))), m(2) * ones(size(XX(:))), XX(:), YY(:)], 'Bandwidth', bandwidth);
%     dens = mvksdensity(results, [XX(:), YY(:), m(3) * ones(size(XX(:))), m(4) * ones(size(XX(:)))], 'Bandwidth', bandwidth);
    dens = reshape(dens, nPoints, nPoints);
    
    % Plots
    contour(XX, YY, dens, 12, 'color', col(i, :))

    i = i + 1;
end
plot(truth(1), truth(2), 'xk', 'markersize', 10)
hL = legend({'$h=0.2$', '$h=0.1$', '$h=0.05$', '$h=0.025$', 'truth'}, 'interpreter', 'laTeX', 'location', 'NE');
box on
yLim = get(gca, 'yLim');
a = get(gca,'XTickLabel');
set(gca, 'XTickLabel', a, 'fontsize', fontsizeTICK)
b = get(gca, 'YTick');
set(gca, 'YTickLabel', cellstr(num2str(b', '%.3f')), 'fontsize', fontsizeTICK)
xlabel('$q_1$', 'interpreter', 'laTeX', 'fontsize', fontsizeLAB)
ylabel('$q_2$', 'interpreter', 'laTeX', 'fontsize', fontsizeLAB)
xLim = get(gca, 'xLim');
yLim = get(gca, 'yLim');

print -depsc2 ../../../Reports/ModelError/BayesDet.eps
close

createFigure(W, H, 'enhanced',enhanced);
hold on
i = 1;
for fileName = {'ODE_IC_ERREST_KEPLER_h_02.txt', 'ODE_IC_ERREST_KEPLER_h_01.txt', 'ODE_IC_ERREST_KEPLER_h_005.txt', 'ODE_IC_ERREST_KEPLER_h_0025.txt', 'ODE_IC_DET_KEPLER_h_small.txt'}

    results = dlmread(fileName{1});
   
    % Evaluate 2D density
    nPoints = 50;
    xEval = linspace(min(results(:, 3)), max(results(:, 3)), nPoints);
    yEval = linspace(min(results(:, 4)), max(results(:, 4)), nPoints);
    [XX, YY] = meshgrid(xEval, yEval);
    
    [n, d] = size(results);    
    bandwidth = std(results) * ((4 / (n * (d + 2)))^(1 / (d + 4)));
    m = mean(results);
        
    dens = mvksdensity(results, [m(1) * ones(size(XX(:))), m(2) * ones(size(XX(:))), XX(:), YY(:)], 'Bandwidth', bandwidth);
%     dens = mvksdensity(results, [XX(:), YY(:), m(3) * ones(size(XX(:))), m(4) * ones(size(XX(:)))], 'Bandwidth', bandwidth);
    dens = reshape(dens, nPoints, nPoints);
    
    % Plots
    contour(XX, YY, dens, 12, 'color', col(i, :))

    i = i + 1;
end
plot(truth(1), truth(2), 'xk', 'markersize', 10)
hL = legend({'$h=0.2$', '$h=0.1$', '$h=0.05$', '$h=0.025$', 'truth'}, 'interpreter', 'laTeX', 'location', 'NE');
box on
set(gca, 'xLim', xLim)
set(gca, 'yLim', yLim)
a = get(gca,'XTickLabel');
set(gca, 'XTickLabel', a, 'fontsize', fontsizeTICK)
b = get(gca, 'YTick');
set(gca, 'YTickLabel', cellstr(num2str(b', '%.3f')), 'fontsize', fontsizeTICK)
xlabel('$q_1$', 'interpreter', 'laTeX', 'fontsize', fontsizeLAB)
ylabel('$q_2$', 'interpreter', 'laTeX', 'fontsize', fontsizeLAB)

print -depsc2 ../../../Reports/ModelError/BayesEE.eps