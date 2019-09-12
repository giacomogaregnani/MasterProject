clc; clear; close all
%%
load('GOOD_RESULTS_DIFF.mat');

W = 6; H = 6;
fig = createFigure(W, H, 'enhanced', 1);
fontsizeLAB = getLatexTextSize('normalsize', 'enhanced', 1);
fontsizeTICK = getLatexTextSize('small', 'enhanced', 1);

plot(x(:, 1), hom(2) * ones(n), 'k');
hold on
plot(x(:, 1), trueVals(2) * ones(n), 'k--');
for j = 3 : 2 : size(x, 2)
    plot(x(:, 1), x(:, j), 'k-o', 'markersize', 5)
end
set(gca, 'ylim', [0, trueVals(2)+0.1])
xlim([x(1, 1), x(end, 1)])
xlabel('$\varepsilon$', 'interpreter', 'latex')
ylabel('$\hat\Sigma_{N, \delta}(\mathbf z^\varepsilon)$', 'interpreter', 'latex')
% set(gca, 'xtick', [])
% set(gca, 'xticklabel', {'10^{-3}', '10^{-2}', '10^{-1}'})

set(gca, 'fontsize', fontsizeTICK);
set(get(gca, 'xlabel'), 'fontsize', fontsizeLAB);
set(get(gca, 'ylabel'), 'fontsize', fontsizeLAB);

box on
axis square
axpos = get(gca, 'Position');
pause(0.1)
export_fig(fig, '../../../Reports/DraftMultiSDE_19/CaltechNotes/Figures/S.eps', '-nocrop', '-painters')

%%
load('GOOD_RESULTS_DIFF_ZETA.mat');

W = 6; H = 6;
fig = createFigure(W, H, 'enhanced', 1);
fontsizeLAB = getLatexTextSize('normalsize', 'enhanced', 1);
fontsizeTICK = getLatexTextSize('small', 'enhanced', 1);

plot(x(:, 1), hom(2) * ones(n), 'k');
hold on
plot(x(:, 1), trueVals(2) * ones(n), 'k--');
for j = 2 : 2 : size(x, 2)
    semilogx(x(:, 1), x(:, j), 'k-o', 'markersize', 5)
end
set(gca, 'ylim', [0, trueVals(2)+0.1])
xlim([x(1, 1), x(end, 1)])
xlabel('$\zeta$', 'interpreter', 'latex')
ylabel('$\hat \Sigma_{N, \delta}(\mathbf z^\varepsilon)$', 'interpreter', 'latex')

set(gca, 'fontsize', fontsizeTICK);
set(get(gca, 'xlabel'), 'fontsize', fontsizeLAB);
set(get(gca, 'ylabel'), 'fontsize', fontsizeLAB);
set(gca, 'xtick', [0, 2, 3, 4]) 
set(gca, 'xticklabel', {'0', '\beta-1', '\beta', '\beta+1'})

box on
axis square
set(gca, 'Position', axpos)
pause(0.1)
export_fig(fig, '../../../Reports/DraftMultiSDE_19/CaltechNotes/Figures/S_zeta.eps', '-nocrop', '-painters')
