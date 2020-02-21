clc; clear; close all
addpath('resultsKL')


%%

s = 0.8;
b = 1;

alpha = dlmread('ResultsAlpha.txt');
hom = dlmread(['ResultsHom_s' num2str(s,'%.1f') '.txt']);
filter = dlmread(['ResultsFilter_s' num2str(s,'%.1f') '_b' num2str(b) '.txt']);
sub = dlmread(['ResultsSub_s' num2str(s,'%.1f') '.txt']);
noth = dlmread(['ResultsNot_s' num2str(s,'%.1f') '.txt']);

aHom = hom(1:end-1);
zetas = filter(:, end);
nParam = length(alpha);

drift = @(alpha, x) driftKL(x, alpha, nParam);
xTest = -1.4:0.01:1.4;

% figure
% plot(xTest, drift(alpha, xTest));
% hold on
% plot(xTest, drift(aHom, xTest));
% title('truth')

% for i = 1 : length(zetas)

W = 4.5; H = 4.5;
fontsizeLAB = getLatexTextSize('normalsize', 'enhanced', 1);
fontsizeTICK = getLatexTextSize('footnotesize', 'enhanced', 1);
enhanced = 1;

fig = createFigure(W, H, 'enhanced', enhanced);
plot(xTest, drift(alpha, xTest), 'k:');
hold on
plot(xTest, drift(aHom, xTest), 'k--');
plot(xTest, drift(filter(1, 1:end-1), xTest), 'k');
legend({'$\alpha \cdot V^\prime(x)$', '$A \cdot V^\prime(x)$', '$\hat A_{\mathrm{filt}} \cdot V^\prime(x)$'}, 'location', 'N', 'interpreter', 'latex')
xlim([xTest(1), xTest(end)])
xlabel('$x$', 'interpreter', 'latex')
set(gca, 'fontsize', fontsizeTICK);
set(get(gca, 'xlabel'), 'fontsize', fontsizeLAB);
set(get(gca, 'ylabel'), 'fontsize', fontsizeLAB);
export_fig(fig, ['~/Desktop/PaperSDE/Figures/KLFilt.png'], '-painters', '-m5')

fig = createFigure(W, H, 'enhanced', enhanced);
plot(xTest, drift(alpha, xTest), 'k:');
hold on
plot(xTest, drift(aHom, xTest), 'k--');
plot(xTest, drift(sub(1, 1:end-1), xTest), 'k');
legend({'$\alpha \cdot V^\prime(x)$', '$A \cdot V^\prime(x)$', '$\hat A_{\mathrm{subs}} \cdot V^\prime(x)$'}, 'location', 'N', 'interpreter', 'latex')
xlim([xTest(1), xTest(end)])
xlabel('$x$', 'interpreter', 'latex')
set(gca, 'fontsize', fontsizeTICK);
set(get(gca, 'xlabel'), 'fontsize', fontsizeLAB);
set(get(gca, 'ylabel'), 'fontsize', fontsizeLAB);
export_fig(fig, ['~/Desktop/PaperSDE/Figures/KLSubs.png'], '-painters', '-m5')


fig = createFigure(W, H, 'enhanced', enhanced);
plot(xTest, drift(alpha, xTest), 'k:');
hold on
plot(xTest, drift(aHom, xTest), 'k--');
plot(xTest, drift(noth(1, :), xTest), 'k');
legend({'$\alpha \cdot V^\prime(x)$', '$A \cdot V^\prime(x)$', '$\hat A \cdot V^\prime(x)$'}, 'location', 'N', 'interpreter', 'latex')
xlim([xTest(1), xTest(end)])
xlabel('$x$', 'interpreter', 'latex')
set(gca, 'fontsize', fontsizeTICK);
set(get(gca, 'xlabel'), 'fontsize', fontsizeLAB);
set(get(gca, 'ylabel'), 'fontsize', fontsizeLAB);
export_fig(fig, ['~/Desktop/PaperSDE/Figures/KLNothing.png'], '-painters', '-m5')
