clc; clear; close all
%%
sol = dlmread('testParamSol.csv');
sol = [sol(1:length(sol)/2)'; sol(length(sol)/2+1:end)'];
plot(sol');

%%
x = dlmread('testParamAvg.txt');

hom = x(1, 1:2);
x = x(2:end, :);
n = size(x(:, 1));
trueVals = [1, 0.7];

%%
W = 6; H = 6;
fig = createFigure(W, H, 'enhanced', 1);
fontsizeLAB = getLatexTextSize('normalsize', 'enhanced', 1);
fontsizeTICK = getLatexTextSize('small', 'enhanced', 1);

plot(x(:, 1), hom(1) * ones(n), 'k');
hold on
plot(x(:, 1), trueVals(1) * ones(n), 'k--');
for j = 2 : 2 : size(x, 2)
    plot(x(:, 1), x(:, j), 'k-o')
end
% set(gca, 'ylim', [0, trueVals(1)+0.2])
xlim([x(1, 1), x(end, 1)])
xlabel('$\log\varepsilon$', 'interpreter', 'latex')
ylabel('$A$', 'interpreter', 'latex')
set(gca, 'xtick', [1e-3, 1e-2, 1e-1]) 
set(gca, 'xticklabel', {'10^{-3}', '10^{-2}', '10^{-1}'})

set(gca, 'fontsize', fontsizeTICK);
set(get(gca, 'xlabel'), 'fontsize', fontsizeLAB);
set(get(gca, 'ylabel'), 'fontsize', fontsizeLAB);

box on
axis square

%%
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