clc; clear; close all
% addpath('resultsBeta_GOOD')
addpath('resultsBeta')

%%
s = .7;

aHom = dlmread(['ResultsHom_s' num2str(s,'%.1f') '.txt']);
aFil = dlmread(['ResultsFilter_s' num2str(s,'%.1f') '.txt']);

W = 4.5; H = 4.5;
enhanced = 4;
fontsizeLAB = getLatexTextSize('small', 'enhanced', 1);
fontsizeTICK = getLatexTextSize('scriptsize', 'enhanced', 1);
fig = createFigure(W, H, 'enhanced', enhanced);

plot([1,10], [aHom, aHom], 'k--')
hold on
plot(1:10, aFil, 'k-o', 'markersize', 3)

xlabel('$\beta$', 'interpreter', 'laTeX')

set(gca, 'fontsize', fontsizeTICK);
set(get(gca, 'xlabel'), 'fontsize', fontsizeLAB);
set(get(gca, 'ylabel'), 'fontsize', fontsizeLAB);
set(get(gca, 'title'), 'fontsize', fontsizeLAB);
ylim([0, .8])
xlim([1, 10])

title(['$\sigma = $ ' num2str(s)], 'interpreter', 'latex')

box on
niceBox(fig)

% export_fig(fig, ['~/Desktop/PaperSDE/Figures/OUBeta_s' num2str(10*s) '.png'], '-nocrop', '-painters', '-m5')