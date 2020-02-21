clc; clear; close all
addpath('resultsCHEB_GOOD')
% addpath('resultsKL')
%%

s = 1.0;
b = 1;

alpha = dlmread('ResultsAlpha.txt');
hom = dlmread(['ResultsHom_s' num2str(s,'%.1f') '.txt']);
filter = dlmread(['ResultsFilter_s' num2str(s,'%.1f') '_b' num2str(b) '.txt']);
sub = dlmread(['ResultsSub_s' num2str(s,'%.1f') '.txt']);
noth = dlmread(['ResultsNot_s' num2str(s,'%.1f') '.txt']);

aHom = hom(1:end-1);
zetas = filter(:, end);
nParam = length(alpha);

drift = @(alpha, x) driftCheb(x, alpha, nParam);
xTest = -1.65:0.01:1.65;

% figure
% plot(xTest, drift(alpha, xTest));
% hold on
% plot(xTest, drift(aHom, xTest));
% title('truth')

% for i = 1 : length(zetas)

W = 4.2; H = 4.2;
fontsizeLAB = getLatexTextSize('footnotesize', 'enhanced', 1);
fontsizeTICK = getLatexTextSize('scriptsize', 'enhanced', 1);
enhanced = 1;

fig = createFigure(W, H, 'enhanced', enhanced);
plot(xTest, drift(alpha, xTest), 'k:');
hold on
plot(xTest, drift(aHom, xTest), 'k--');
plot(xTest, drift(filter(1, 1:end-1), xTest), 'k');
legend({'$\alpha \cdot V^\prime(x)$', '$A \cdot V^\prime(x)$', '$\widehat A_k^{\varepsilon} \cdot V^\prime(x)$'}, 'location', 'N', 'interpreter', 'latex')
xlim([xTest(1), xTest(end)])
ylim([-3, 25])
xlabel('$x$', 'interpreter', 'latex')
title('Filtering', 'interpreter', 'latex')
set(gca, 'fontsize', fontsizeTICK);
set(get(gca, 'xlabel'), 'fontsize', fontsizeLAB);
set(get(gca, 'ylabel'), 'fontsize', fontsizeLAB);
set(get(gca, 'legend'), 'fontsize', fontsizeLAB);
set(get(gca, 'title'), 'fontsize', fontsizeLAB);

export_fig(fig, ['~/Desktop/PaperSDE/Figures/KLFilt.png'], '-painters', '-m5')

fig = createFigure(W, H, 'enhanced', enhanced);
plot(xTest, drift(alpha, xTest), 'k:');
hold on
plot(xTest, drift(aHom, xTest), 'k--');
plot(xTest, drift(sub(1, 1:end-1), xTest), 'k');
legend({'$\alpha \cdot V^\prime(x)$', '$A \cdot V^\prime(x)$', '$\widehat A_\delta^{\varepsilon} \cdot V^\prime(x)$'}, 'location', 'N', 'interpreter', 'latex')
xlim([xTest(1), xTest(end)])
ylim([-3, 25])
xlabel('$x$', 'interpreter', 'latex')
title('Subsampling', 'interpreter', 'latex')
set(gca, 'fontsize', fontsizeTICK);
set(get(gca, 'xlabel'), 'fontsize', fontsizeLAB);
set(get(gca, 'ylabel'), 'fontsize', fontsizeLAB);
set(get(gca, 'legend'), 'fontsize', fontsizeLAB);
set(get(gca, 'title'), 'fontsize', fontsizeLAB);

export_fig(fig, ['~/Desktop/PaperSDE/Figures/KLSubs.png'], '-painters', '-m5')


fig = createFigure(W, H, 'enhanced', enhanced);
plot(xTest, drift(alpha, xTest), 'k:');
hold on
plot(xTest, drift(aHom, xTest), 'k--');
plot(xTest, drift(noth(1, :), xTest), 'k');
legend({'$\alpha \cdot V^\prime(x)$', '$A \cdot V^\prime(x)$', '$\widehat A^{\varepsilon} \cdot V^\prime(x)$'}, 'location', 'N', 'interpreter', 'latex')
xlim([xTest(1), xTest(end)])
ylim([-3, 25])
xlabel('$x$', 'interpreter', 'latex')
title('No preprocessing', 'interpreter', 'latex')
set(gca, 'fontsize', fontsizeTICK);
set(get(gca, 'xlabel'), 'fontsize', fontsizeLAB);
set(get(gca, 'ylabel'), 'fontsize', fontsizeLAB);
set(get(gca, 'legend'), 'fontsize', fontsizeLAB);
set(get(gca, 'title'), 'fontsize', fontsizeLAB);

export_fig(fig, ['~/Desktop/PaperSDE/Figures/KLNothing.png'], '-painters', '-m5')