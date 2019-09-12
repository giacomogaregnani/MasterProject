clc; clear; close all
%%
load('GOOD_RESULTS_DRIFT.mat');

W = 6; H = 6;
fig = createFigure(W, H, 'enhanced', 1);
fontsizeLAB = getLatexTextSize('normalsize', 'enhanced', 1);
fontsizeTICK = getLatexTextSize('small', 'enhanced', 1);

semilogx(x(:, 1), hom(1) * ones(n), 'k');
hold on
semilogx(x(:, 1), trueVals(1) * ones(n), 'k--');
for j = 2 : 2 : size(x, 2)
    semilogx(x(:, 1), x(:, j), 'k-o', 'markersize', 5)
end
set(gca, 'ylim', [0, trueVals(1)+0.2])
xlim([x(1, 1), x(end, 1)])
xlabel('$\log\varepsilon$', 'interpreter', 'latex')
ylabel('$\hat A_{N, \delta}(\mathbf z^\varepsilon)$', 'interpreter', 'latex')
set(gca, 'xtick', [1e-3, 1e-2, 1e-1]) 
set(gca, 'xticklabel', {'10^{-3}', '10^{-2}', '10^{-1}'})
set(gca, 'ytick', [0, 0.2, 0.4, 0.6, 0.8, 1, 1.2])

set(gca, 'fontsize', fontsizeTICK);
set(get(gca, 'xlabel'), 'fontsize', fontsizeLAB);
set(get(gca, 'ylabel'), 'fontsize', fontsizeLAB);

box on
axis square
axpos = get(gca, 'Position');
pause(0.1)
export_fig(fig, '../../../Reports/DraftMultiSDE_19/CaltechNotes/Figures/A.eps', '-nocrop', '-painters')

%%
load('GOOD_RESULTS_DRIFT_ZETA.mat');

W = 6; H = 6;
fig = createFigure(W, H, 'enhanced', 1);
fontsizeLAB = getLatexTextSize('normalsize', 'enhanced', 1);
fontsizeTICK = getLatexTextSize('small', 'enhanced', 1);

plot(x(:, 1), hom(1) * ones(n), 'k');
hold on
plot(x(:, 1), trueVals(1) * ones(n), 'k--');
for j = 2 : 2 : size(x, 2)
    semilogx(x(:, 1), x(:, j), 'k-o', 'markersize', 5)
end
set(gca, 'ylim', [0, trueVals(1)+0.2])
xlim([x(end, 1), x(1, 1)])
xlabel('$\zeta$', 'interpreter', 'latex')
ylabel('$\hat A_{N, \delta}(\mathbf z^\varepsilon)$', 'interpreter', 'latex')

set(gca, 'xtick', [2, 2.5, 3])
set(gca, 'xticklabel', {'\beta-1', '\beta-1/2', '\beta'});
set(gca, 'ytick', [0, 0.2, 0.4, 0.6, 0.8, 1, 1.2])

set(gca, 'fontsize', fontsizeTICK);
set(get(gca, 'xlabel'), 'fontsize', fontsizeLAB);
set(get(gca, 'ylabel'), 'fontsize', fontsizeLAB);

box on
axis square
set(gca, 'Position', axpos)
pause(0.1)
export_fig(fig, '../../../Reports/DraftMultiSDE_19/CaltechNotes/Figures/A_zeta.eps', '-nocrop', '-painters')
