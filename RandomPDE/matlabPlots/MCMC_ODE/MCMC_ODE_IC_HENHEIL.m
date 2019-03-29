clc; clear; close all
addpath('RESULTSHENONSTORMVERL');
%%

load('axpos')

truth = [0, 0.1];
col = colormap(parula(8));
close
col = col(1:2:end, :);
enhanced = 1;
fontsizeLAB = getLatexTextSize('normalsize', 'enhanced', enhanced);
fontsizeTICK = getLatexTextSize('footnotesize', 'enhanced', enhanced);
W = 4; H = 4;

createFigure(W, H, 'enhanced',enhanced);
hold on
i = 1;
for fileName = {'ODE_IC_PROB_KEPLER_h_02.txt', 'ODE_IC_PROB_KEPLER_h_01.txt', 'ODE_IC_PROB_KEPLER_h_005.txt', 'ODE_IC_PROB_KEPLER_h_0025.txt'}

    results = dlmread(fileName{1});
   
    % Evaluate 2D density
    nPoints = 50;
    xEval = linspace(min(results(:, 1)), max(results(:, 1)), nPoints);
    yEval = linspace(min(results(:, 2)), max(results(:, 2)), nPoints);
    [XX, YY] = meshgrid(xEval, yEval);
    
    [n, d] = size(results);
    bandwidth = std(results) * ((4 / (n * (d + 2)))^(1 / (d + 4)));
    m = mean(results);
        
    dens = mvksdensity(results, [XX(:), YY(:), m(3) * ones(size(XX(:))), m(4) * ones(size(XX(:)))], 'Bandwidth', bandwidth);
    dens = reshape(dens, nPoints, nPoints);
    
    % Plots
    contour(XX, YY, dens, 12, 'color', col(i, :))

    i = i + 1;
end
plot(truth(1), truth(2), 'xk', 'markersize', 10)

box on
yLim = [0.09, 0.12];
yTick = [0.09, 0.1, 0.11, 0.12];
set(gca, 'yLim', yLim);
set(gca, 'yTick', yTick);
xLim = get(gca, 'xLim');

xlabel('$v_1$', 'interpreter', 'laTeX', 'fontsize', fontsizeLAB)
ylabel('$v_2$', 'interpreter', 'laTeX', 'fontsize', fontsizeLAB)
set(gca, 'innerposition', axpos)

print -depsc2 ../../../Reports/PaperRTSRK_18/VERSION12/BayesProb.eps
close

% figure
createFigure(W, H, 'enhanced',enhanced);
hold on
i = 1;
for fileName = {'ODE_IC_DET_EXP_KEPLER_h_02.txt', 'ODE_IC_DET_EXP_KEPLER_h_01.txt', 'ODE_IC_DET_EXP_KEPLER_h_005.txt', 'ODE_IC_DET_EXP_KEPLER_h_0025.txt'}

    results = dlmread(fileName{1});
   
    % Evaluate 2D density
    nPoints = 50;
    xEval = linspace(min(results(:, 1)), max(results(:, 1)), nPoints);
    yEval = linspace(min(results(:, 2)), max(results(:, 2)), nPoints);
    [XX, YY] = meshgrid(xEval, yEval);
    
    [n, d] = size(results);    
    bandwidth = std(results) * ((4 / (n * (d + 2)))^(1 / (d + 4)));
    m = mean(results);
        
    dens = mvksdensity(results, [XX(:), YY(:), m(3) * ones(size(XX(:))), m(4) * ones(size(XX(:)))], 'Bandwidth', bandwidth);
    dens = reshape(dens, nPoints, nPoints);
    
    % Plots
    contour(XX, YY, dens, 12, 'color', col(i, :))

    i = i + 1;
end
plot(truth(1), truth(2), 'xk', 'markersize', 10)

box on
set(gca, 'xLim', xLim);
set(gca, 'yLim', yLim);
xlabel('$v_1$', 'interpreter', 'laTeX', 'fontsize', fontsizeLAB)
ylabel('$v_2$', 'interpreter', 'laTeX', 'fontsize', fontsizeLAB)
set(gca, 'innerposition', axpos)

print -depsc2 ../../../Reports/PaperRTSRK_18/VERSION12/BayesHeun.eps
close

createFigure(W, H, 'enhanced',enhanced);
hold on
i = 1;
for fileName = {'ODE_IC_DET_KEPLER_h_02.txt', 'ODE_IC_DET_KEPLER_h_01.txt', 'ODE_IC_DET_KEPLER_h_005.txt', 'ODE_IC_DET_KEPLER_h_0025.txt'}
    
    results = dlmread(fileName{1});
   
%     Evaluate 2D density
    nPoints = 50;
    xEval = linspace(min(results(:, 1)), max(results(:, 1)), nPoints);
    yEval = linspace(min(results(:, 2)), max(results(:, 2)), nPoints);
    [XX, YY] = meshgrid(xEval, yEval);
    
    [n, d] = size(results);
    bandwidth = std(results) * ((4 / (n * (d + 2)))^(1 / (d + 4)));
    m = mean(results);
        
    dens = mvksdensity(results, [XX(:), YY(:), m(3) * ones(size(XX(:))), m(4) * ones(size(XX(:)))], 'Bandwidth', bandwidth);
    dens = reshape(dens, nPoints, nPoints);
    
%     Plots
    contour(XX, YY, dens, 12, 'color', col(i, :))

    i = i + 1;
end
plot(truth(1), truth(2), 'xk', 'markersize', 10)

box on
set(gca, 'xLim', xLim);
set(gca, 'yLim', yLim);
xlabel('$v_1$', 'interpreter', 'laTeX', 'fontsize', fontsizeLAB)
ylabel('$v_2$', 'interpreter', 'laTeX', 'fontsize', fontsizeLAB)
set(gca, 'innerposition', axpos)

print -depsc2 ../../../Reports/PaperRTSRK_18/VERSION12/BayesDet.eps

