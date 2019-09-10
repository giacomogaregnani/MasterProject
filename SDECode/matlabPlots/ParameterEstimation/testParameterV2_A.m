clc; clear; close all
%%
x = dlmread('testParamZeta.txt');

hom = x(1, 1:2);
x = x(2:end, :);
n = size(x(:, 1));
trueVals = [1, 0.5];

%%
W = 6; H = 6;
fig = createFigure(W, H, 'enhanced', 1);
fontsizeLAB = getLatexTextSize('normalsize', 'enhanced', 1);
fontsizeTICK = getLatexTextSize('small', 'enhanced', 1);

plot(x(:, 1), hom(1) * ones(n), 'k');
hold on
plot(x(:, 1), trueVals(1) * ones(n), 'k--');
for j = 2 : 2 : size(x, 2)
    semilogx(x(:, 1), x(:, j), 'k-o')
end
set(gca, 'ylim', [0, trueVals(1)+0.2])
xlim([x(end, 1), x(1, 1)])
xlabel('$\zeta$', 'interpreter', 'latex')
ylabel('$A$', 'interpreter', 'latex')
set(gca, 'xtick', [2, 2.5, 3])
set(gca, 'xticklabel', {'\beta-1', '\beta-1/2', '\beta'});

set(gca, 'fontsize', fontsizeTICK);
set(get(gca, 'xlabel'), 'fontsize', fontsizeLAB);
set(get(gca, 'ylabel'), 'fontsize', fontsizeLAB);

box on
axis square
export_fig(fig, '../../../Reports/DraftMultiSDE_19/CaltechNotes/Figures/A_zeta.eps', '-nocrop', '-painters')