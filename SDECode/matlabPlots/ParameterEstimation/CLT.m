clc; clear; close all
addpath('resultsCLT');
%%

T = 50;

aHom = dlmread(['OUHom', num2str(T), '.txt']);
A = dlmread(['OU', num2str(T), '.txt']);

W = 13; H = 4.5;
enhanced = 3;
fontsizeLAB = getLatexTextSize('small', 'enhanced', 1);
fontsizeTICK = getLatexTextSize('scriptsize', 'enhanced', 1);
% fig = figure;
fig = createFigure(W, H, 'enhanced', enhanced);

alpha = 1;
vv = 2 * alpha * (1 + aHom);

mEmp = mean(sqrt(T) * (A - aHom));
sEmp = std(sqrt(T) * A);
mTh = 0;
% sTh = sqrt(vv);
sTh = sqrt(dlmread(['OUVar', num2str(T), '.txt']));
histogram(sqrt(T) * (A - aHom), 'normalization', 'pdf', 'numbins', 30, 'faceColor', [0.8, 0.8, 0.8])
hold on
x = linspace(norminv(0.001, mEmp, sEmp), norminv(0.999, mEmp, sEmp), 100);
plot(x, normpdf(x, mEmp, sEmp), 'k')
x = linspace(norminv(0.001, mTh, sTh), norminv(0.999, mTh, sTh), 100);
plot(x, normpdf(x, mTh, sTh), 'k--')

xLim = get(gca, 'xlim');
maxxLim = max(xLim);
set(gca, 'xlim', [-maxxLim, maxxLim])

set(gca, 'fontsize', fontsizeTICK);
set(gca, 'ytick', [])

legend({'$\sqrt{T}(\widehat A_k^\varepsilon(T) - A)$', '$\widehat \rho_A$', '$\rho_A$'}, 'interpreter', 'latex', 'fontsize', fontsizeLAB);
legend({'$\sqrt{T}(\widehat A_k^\varepsilon(T) - A)$', 'Gaussian fit', 'Theory'}, 'interpreter', 'latex', 'fontsize', fontsizeLAB);

placeLegendExtreme(fig, 3, 'NW');
niceBox(fig)
          
% export_fig(fig, ['~/Desktop/PaperSDE/Figures/CLT.png'], '-nocrop', '-painters', '-m5')