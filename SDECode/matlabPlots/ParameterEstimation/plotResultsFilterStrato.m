clc; clear; close all
%%

aFil = dlmread('ResultsFilterStrato.txt');
hom = dlmread('ResultsFilterHomStrato.txt');
aFul = dlmread('ResultsFilterFullStrato.txt');
aSub = dlmread('ResultsFilterSubsStrato.txt');
aHom = hom(1);

W = 6; H = 6;
fig = createFigure(W, H, 'enhanced', 1);
fontsizeLAB = getLatexTextSize('normalsize', 'enhanced', 1);
fontsizeTICK = getLatexTextSize('small', 'enhanced', 1);

f = zeros(3, 1);
% l = Line.array()
[h(1), l(1)] = histAndPDF(aFil, 'FaceColor', 'k', 'NumBins', 20, 'FaceAlpha', 0.25, 'LineStyle', '--', 'LineColor', 'k', 'LineWidth', .5);
[h(2), l(2)] = histAndPDF(aSub, 'FaceColor', 'r', 'NumBins', 20, 'FaceAlpha', 0.25, 'LineColor', 'k', 'LineWidth', .5);
% [h(3), l(3)] = histAndPDF(aFul, 'FaceColor', 'b', 'NumBins', 20, 'FaceAlpha', 0.25, '--', 'LineColor', 'k', 'LineWidth', .5);
% legend('Filter', 'Full', 'Subsampling', 'Location', 'best')

lMax = max([l(:).YData]);
plot([aHom, aHom], [0, lMax], 'k--', 'LineWidth', .5)
legend([l(1),l(2)], {'filter', 'subsamp'}, 'Location', 'NE')

xlabel('$A$', 'interpreter', 'laTeX')

set(gca, 'fontsize', fontsizeTICK);
set(get(gca, 'xlabel'), 'fontsize', fontsizeLAB);
set(get(gca, 'ylabel'), 'fontsize', fontsizeLAB);

% xLim = get(gca, 'xLim');
% set(gca, 'xLim', [aHom-0.05, xLim(2)])
set(gca, 'yTick', [])
set(gca, 'yTickLabel', [])
box on
axis square
export_fig(fig, '../../../Reports/DraftMultiSDE_19/CaltechNotes_2/Figures/A_Strato_OU_S1.png', '-nocrop', '-painters', '-m10')