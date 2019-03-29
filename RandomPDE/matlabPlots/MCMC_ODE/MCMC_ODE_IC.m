clc; clear; close all
%%

truth = [1.5, -pi];
col = 'rbg';
enhanced = 1;
fontsizeLAB = getLatexTextSize('normalsize', 'enhanced', enhanced);
fontsizeTICK = getLatexTextSize('footnotesize', 'enhanced', enhanced);
W = 7; H = 7;

createFigure(W, H, 'enhanced',enhanced);
hold on
i = 1;
for fileName = {'ODE_IC_PROB_PEND_h_01.txt', 'ODE_IC_PROB_PEND_h_005.txt', 'ODE_IC_PROB_PEND_h_0025.txt'}
    
    results = dlmread(fileName{1});
   
    % Evaluate 2D density
    nPoints = 60;
    xEval = linspace(0.995*min(results(:, 1)), 1.005*max(results(:, 1)), nPoints);
    yEval = linspace(1.005*min(results(:, 2)), 0.995*max(results(:, 2)), nPoints);
    [XX, YY] = meshgrid(xEval, yEval);
    
%     dens = akde([results(:, 1), results(:, 2)], [XX(:), YY(:)]);
    dens = mvksdensity([results(:, 1), results(:, 2)], [XX(:), YY(:)]);
    dens = reshape(dens, nPoints, nPoints);
    
    % Plots
    contour(XX, YY, dens, 10, 'color', col(i))

    i = i + 1;
end
plot(truth(1), truth(2), 'xk', 'markersize', 10)
hL = legend({'$h=0.1$', '$h=0.05$', '$h=0.025$', 'truth'}, 'interpreter', 'laTeX', 'location', 'SW');
currentLegendPosition = hL.Position;
newLegendPosition = [0.18 0.115 currentLegendPosition([3 4])];
hL.Position = newLegendPosition;
box on
xLim = get(gca, 'xLim');
yLim = get(gca, 'yLim');

a = get(gca,'XTickLabel');
set(gca, 'XTickLabel', a, 'fontsize', fontsizeTICK)
a = get(gca,'YTickLabel');
set(gca, 'YTickLabel', a, 'fontsize', fontsizeTICK)
xlabel('$p_0$', 'interpreter', 'laTeX', 'fontsize', fontsizeLAB)
ylabel('$q_0$', 'interpreter', 'laTeX', 'fontsize', fontsizeLAB)
title('Probabilistic solver', 'interpreter', 'latex', 'fontsize', fontsizeLAB)
% print -depsc2 ../../../Reports/RTSRK_Hamiltonian/Hamiltonian3/BayesProb.eps


close

createFigure(W, H, 'enhanced',enhanced);
hold on
i = 1;
for fileName = {'ODE_IC_DET_PEND_h_01.txt', 'ODE_IC_DET_PEND_h_005.txt', 'ODE_IC_DET_PEND_h_0025.txt'}
    
    results = dlmread(fileName{1});
   
    % Evaluate 2D density
    nPoints = 60;
    xEval = linspace(min(results(:, 1)), max(results(:, 1)), nPoints);
    yEval = linspace(min(results(:, 2)), max(results(:, 2)), nPoints);
    [XX, YY] = meshgrid(xEval, yEval);
    
    dens = mvksdensity([results(:, 1), results(:, 2)], [XX(:), YY(:)]);
    dens = reshape(dens, nPoints, nPoints);
    
    % Plots
    contour(XX, YY, dens, 10, 'color', col(i))

    i = i + 1;
end
plot(truth(1), truth(2), 'xk', 'markersize', 10)
hL = legend({'$h=0.1$', '$h=0.05$', '$h=0.025$', 'truth'}, 'interpreter', 'laTeX', 'location', 'NW');
currentLegendPosition = hL.Position;
newLegendPosition = [0.18 0.115 currentLegendPosition([3 4])];
hL.Position = newLegendPosition;
box on
set(gca, 'xLim', xLim);
set(gca, 'yLim', yLim);
a = get(gca,'XTickLabel');
set(gca, 'XTickLabel', a, 'fontsize', fontsizeTICK)
a = get(gca,'YTickLabel');
set(gca, 'YTickLabel', a, 'fontsize', fontsizeTICK)
xlabel('$p_0$', 'interpreter', 'laTeX', 'fontsize', fontsizeLAB)
ylabel('$q_0$', 'interpreter', 'laTeX', 'fontsize', fontsizeLAB)
title('Deterministic solver', 'interpreter', 'latex', 'fontsize', fontsizeLAB)
% print -depsc2 ../../../Reports/RTSRK_Hamiltonian/Hamiltonian3/BayesDet.eps