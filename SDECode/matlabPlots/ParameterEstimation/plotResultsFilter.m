clc; clear; close all
%%

aFil = dlmread('ResultsFilter.txt');
hom = dlmread('ResultsFilterHom.txt');
aFul = dlmread('ResultsFilterFull.txt');
aSub = dlmread('ResultsFilterSubs.txt');
aHom = hom(1);

W = 6; H = 6;
fig = createFigure(W, H, 'enhanced', 1);
fontsizeLAB = getLatexTextSize('normalsize', 'enhanced', 1);
fontsizeTICK = getLatexTextSize('small', 'enhanced', 1);

f = zeros(3, 1);
% l = Line.array()
[h(1), l(1)] = histAndPDF(aFil, 'FaceColor', 'k', 'NumBins', 20, 'FaceAlpha', 0.25, '--', 'LineColor', 'k', 'LineWidth', .5);
[h(2), l(2)] = histAndPDF(aSub, 'FaceColor', 'r', 'NumBins', 20, 'FaceAlpha', 0.25, '--', 'LineColor', 'k', 'LineWidth', .5);
% [h(3), l(3)] = histAndPDF(aFul, 'FaceColor', 'b', 'NumBins', 20, 'FaceAlpha', 0.25, '--', 'LineColor', 'k', 'LineWidth', .5);
% legend('Filter', 'Full', 'Subsampling', 'Location', 'best')

l = max([l(:).YData]);
plot([aHom, aHom], [0, l], 'k--', 'LineWidth', .5)
% legend(h, {'fil', 'ful', 'sub'})

set(gca, 'fontsize', fontsizeTICK);
set(get(gca, 'xlabel'), 'fontsize', fontsizeLAB);
set(get(gca, 'ylabel'), 'fontsize', fontsizeLAB);

xLim = get(gca, 'xLim');
% set(gca, 'xLim', [aHom-0.05, xLim(2)])
set(gca, 'yTick', [])
set(gca, 'yTickLabel', [])
box on
axis square
% print -depsc '../../../Reports/DraftMultiSDE_19/CaltechNotes_2/Figures/A_OU.eps'
%export_fig(fig, '../../../Reports/DraftMultiSDE_19/CaltechNotes_2/Figures/A_OU.png', '-nocrop', '-painters', '-m10')

%% Diffusion

sFil = dlmread('ResultsFilterDiff.txt');
sFul = dlmread('ResultsFilterFullDiff.txt');
sSub = dlmread('ResultsFilterSubsDiff.txt');
sHom = hom(2);

W = 6; H = 6;
fig = createFigure(W, H, 'enhanced', 1);
fontsizeLAB = getLatexTextSize('normalsize', 'enhanced', 1);
fontsizeTICK = getLatexTextSize('small', 'enhanced', 1);

f = zeros(3, 1);
clear l
[h(1), l(1)] = histAndPDF(sFil, 'FaceColor', 'k', 'NumBins', 20, 'FaceAlpha', 0.25, '--', 'LineColor', 'k', 'LineWidth', .5);
[h(2), l(2)] = histAndPDF(sSub, 'FaceColor', 'r', 'NumBins', 20, 'FaceAlpha', 0.25, '--', 'LineColor', 'k', 'LineWidth', .5);
% [h(3), l(3)] = histAndPDF(sFul, 'FaceColor', 'b', 'NumBins', 20, 'FaceAlpha', 0.25, '--', 'LineColor', 'k', 'LineWidth', .5);
% legend('Filter', 'Full', 'Subsampling', 'Location', 'best')

l = max([l(:).YData]);
plot([sHom, sHom], [0, l], 'k--', 'LineWidth', .5)
% legend(h, {'fil', 'ful', 'sub'})

set(gca, 'fontsize', fontsizeTICK);
set(get(gca, 'xlabel'), 'fontsize', fontsizeLAB);
set(get(gca, 'ylabel'), 'fontsize', fontsizeLAB);

xLim = get(gca, 'xLim');
% set(gca, 'xLim', [aHom-0.05, xLim(2)])
set(gca, 'yTick', [])
set(gca, 'yTickLabel', [])
box on
axis square
% print -depsc '../../../Reports/DraftMultiSDE_19/CaltechNotes_2/Figures/A_OU.eps'
%export_fig(fig, '../../../Reports/DraftMultiSDE_19/CaltechNotes_2/Figures/S_OU.png', '-nocrop', '-painters', '-m10')