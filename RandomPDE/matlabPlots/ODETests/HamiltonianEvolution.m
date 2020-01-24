clc; clear; close all
addpath('resultsHamiltonianTraj')

%%
results = dlmread('test_h_02.txt');
meanResult = mean(results(:, 3:end), 2);
stdResult = std(results(:, 3:end), 0, 2);

fontsizeLAB = getLatexTextSize('normalsize', 'enhanced', 1);
fontsizeTICK = getLatexTextSize('footnotesize', 'enhanced', 1);
W = 12; H = 5;

createFigure(W, H, 'enhanced',1);
hold on

plot(results(:, 1), meanResult + 2.0 * stdResult, 'color', 0.7 * ones(1, 3))
plot(results(:, 1), meanResult - 2.0 * stdResult, 'color' , 0.7 * ones(1, 3))
plot(results(:, 1), meanResult, 'k')
plot([results(1, 1), results(end, 1)], [results(1, 2), results(1, 2)], 'k--')

box on

a = get(gca,'XTickLabel');
set(gca, 'XTickLabel', a, 'fontsize', fontsizeTICK)
a = get(gca,'YTickLabel');
set(gca, 'YTickLabel', a, 'fontsize', fontsizeTICK)

xlabel('$t_n$', 'interpreter', 'latex', 'fontsize', fontsizeLAB)
ylabel('$Q(Y_n)$', 'interpreter', 'latex', 'fontsize', fontsizeLAB)

yLim = get(gca, 'yLim');

print -depsc2 ../../../Reports/PaperRTSRK_18/VERSION17/Paper/HamiltonianTraj.eps

%%

results = dlmread('test_h_01.txt');
meanResult = mean(results(:, 3:end), 2);
stdResult = std(results(:, 3:end), 0, 2);

createFigure(W, H, 'enhanced',1);
hold on

plot(results(:, 1), meanResult + 2.0 * stdResult, 'color', 0.7 * ones(1, 3))
plot(results(:, 1), meanResult - 2.0 * stdResult, 'color' , 0.7 * ones(1, 3))
plot(results(:, 1), meanResult, 'k')
plot([results(1, 1), results(end, 1)], [results(1, 2), results(1, 2)], 'k--')

box on

set(gca, 'yLim', yLim)

a = get(gca,'XTickLabel');
set(gca, 'XTickLabel', a, 'fontsize', fontsizeTICK)
a = get(gca,'YTickLabel');
set(gca, 'YTickLabel', a, 'fontsize', fontsizeTICK)

xlabel('$t_n$', 'interpreter', 'latex', 'fontsize', fontsizeLAB)
ylabel('$Q(Y_n)$', 'interpreter', 'latex', 'fontsize', fontsizeLAB)


print -depsc2 ../../../Reports/PaperRTSRK_18/VERSION17/Paper/HamiltonianTraj2.eps
