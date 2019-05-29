clc; clear; % close all
addpath('RESULTSHENONSTORMVERL');
%%

truth = [0.5, 0]; 
col = colormap(parula(8));
% close
col = col(1:2:end, :);
enhanced = 1;
fontsizeLAB = getLatexTextSize('normalsize', 'enhanced', enhanced);
fontsizeTICK = getLatexTextSize('footnotesize', 'enhanced', enhanced);
W = 4; H = 4;

createFigure(W, H, 'enhanced',enhanced);
hold on
i = 1;
for fileName = {'ODE_IC_ADD_KEPLER_REVIEW_h_02.txt', 'ODE_IC_ADD_KEPLER_REVIEW_h_01.txt', 'ODE_IC_ADD_KEPLER_REVIEW_h_005.txt', 'ODE_IC_ADD_KEPLER_REVIEW_h_0025.txt'}

    results = dlmread(fileName{1});
   
    % Evaluate 2D density
    nPoints = 20;
    xEval = linspace(min(results(:, 3)), max(results(:, 3)), nPoints);
    yEval = linspace(min(results(:, 4)), max(results(:, 4)), nPoints);
    [XX, YY] = meshgrid(xEval, yEval);
    
    [n, d] = size(results);    
    bandwidth = std(results) * ((4 / (n * (d + 2)))^(1 / (d + 4)));
    m = mean(results);
        
    dens = mvksdensity(results, [m(1) * ones(size(XX(:))), m(2) * ones(size(XX(:))), XX(:), YY(:)], 'Bandwidth', bandwidth);
    dens = reshape(dens, nPoints, nPoints);
    
    % Plots
    contour(XX, YY, dens, 12, 'color', col(i, :))

    i = i + 1;
end
plot(truth(1), truth(2), 'xk', 'markersize', 10)

box on
xLim = [0.485, 0.51]; 
yLim = [-8, 4] * 1e-3;
set(gca, 'xLim', xLim);
set(gca, 'yLim', yLim);
xlabel('$w_1$', 'interpreter', 'laTeX', 'fontsize', fontsizeLAB)
ylabel('$w_2$', 'interpreter', 'laTeX', 'fontsize', fontsizeLAB)
scale = 0.08;
axpos = get(gca, 'Position');
axpos(2) = axpos(2)+scale*axpos(4);
axpos(4) = (1-scale)*axpos(4);
set(gca, 'innerposition', axpos);

%print -depsc2 ../../../Reports/PaperRTSRK_18/VERSION12/BayesHeun2.eps
% close

createFigure(W, H, 'enhanced',enhanced);
hold on
i = 1;
for fileName = {'ODE_IC_PROB_KEPLER_REVIEW_h_02.txt', 'ODE_IC_PROB_KEPLER_REVIEW_h_01.txt', 'ODE_IC_PROB_KEPLER_REVIEW_h_005.txt', 'ODE_IC_PROB_KEPLER_REVIEW_h_0025.txt'}

    results = dlmread(fileName{1});
   
    % Evaluate 2D density
    nPoints = 20;
    xEval = linspace(min(results(:, 3)), max(results(:, 3)), nPoints);
    yEval = linspace(min(results(:, 4)), max(results(:, 4)), nPoints);
    [XX, YY] = meshgrid(xEval, yEval);
    
    [n, d] = size(results);
    bandwidth = std(results) * ((4 / (n * (d + 2)))^(1 / (d + 4)));
    m = mean(results);
        
    dens = mvksdensity(results, [m(1) * ones(size(XX(:))), m(2) * ones(size(XX(:))), XX(:), YY(:)], 'Bandwidth', bandwidth);    
    dens = reshape(dens, nPoints, nPoints);
    
    % Plots
    contour(XX, YY, dens, 12, 'color', col(i, :))

    i = i + 1;
end
plot(truth(1), truth(2), 'xk', 'markersize', 10)

box on
set(gca, 'xLim', xLim);
set(gca, 'yLim', yLim);
xlabel('$w_1$', 'interpreter', 'laTeX', 'fontsize', fontsizeLAB)
ylabel('$w_2$', 'interpreter', 'laTeX', 'fontsize', fontsizeLAB)
set(gca, 'innerposition', axpos)

%print -depsc2 ../../../Reports/PaperRTSRK_18/VERSION12/BayesProb2.eps
% close

createFigure(W, H, 'enhanced',enhanced);
hold on
i = 1;
for fileName = {'ODE_IC_DET_KEPLER_REVIEW_h_02.txt', 'ODE_IC_DET_KEPLER_REVIEW_h_01.txt', 'ODE_IC_DET_KEPLER_REVIEW_h_005.txt', 'ODE_IC_DET_KEPLER_REVIEW_h_0025.txt'}
    
    results = dlmread(fileName{1});
   
%     Evaluate 2D density
    nPoints = 20;
    xEval = linspace(min(results(:, 3)), max(results(:, 3)), nPoints);
    yEval = linspace(min(results(:, 4)), max(results(:, 4)), nPoints);
    [XX, YY] = meshgrid(xEval, yEval);
    
    [n, d] = size(results);
    bandwidth = std(results) * ((4 / (n * (d + 2)))^(1 / (d + 4)));
    m = mean(results);
        
    dens = mvksdensity(results, [m(1) * ones(size(XX(:))), m(2) * ones(size(XX(:))), XX(:), YY(:)], 'Bandwidth', bandwidth);    
    dens = reshape(dens, nPoints, nPoints);
    
%     Plots
    contour(XX, YY, dens, 12, 'color', col(i, :))

    i = i + 1;
end
plot(truth(1), truth(2), 'xk', 'markersize', 10)

box on
set(gca, 'xLim', xLim);
set(gca, 'yLim', yLim);
xlabel('$w_1$', 'interpreter', 'laTeX', 'fontsize', fontsizeLAB)
ylabel('$w_2$', 'interpreter', 'laTeX', 'fontsize', fontsizeLAB)
set(gca, 'innerposition', axpos)

%print -depsc2 ../../../Reports/PaperRTSRK_18/VERSION12/BayesDet2.eps

save('axpos', 'axpos')