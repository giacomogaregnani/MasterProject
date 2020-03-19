clc; clear; close all
% addpath('resultsOU_epsSmall')
% addpath('resultsOU_GOOD')
% addpath('resultsQU_epsBig2')
addpath('resultsOU_OtherFilter')

%%
s = 1;
zetas = dlmread('zetas.txt');

W = 4.5; H = 4.5;
% W = 3.5; H = 3.5;
enhanced = 1;

aHom = dlmread(['ResultsHom_s' num2str(s,'%.1f') '.txt']);
aSub = dlmread(['ResultsSub_s' num2str(s,'%.1f') '.txt']);
aSub_m = aSub(1:ceil(length(aSub)/2));
aSub_s = aSub(ceil(length(aSub)/2)+1:end);

fig = createFigure(W, H, 'enhanced', enhanced);
fontsizeLAB = getLatexTextSize('small', 'enhanced', 1);
fontsizeTICK = getLatexTextSize('scriptsize', 'enhanced', 1);

plot([0, 1], [aHom, aHom], 'k--')
hold on
plot(zetas, aSub_m, 'k-o', 'markersize', 3)
% plot(zetas, aSub_m + 2*aSub_s, 'k--')
% plot(zetas, aSub_m - 2*aSub_s, 'k--')

xlabel('$\zeta$', 'interpreter', 'laTeX')

title('Subsampling', 'interpreter', 'latex')
% title(['$\sigma = ', num2str(s), '$'], 'interpreter', 'latex')

set(gca, 'fontsize', fontsizeTICK);
set(get(gca, 'xlabel'), 'fontsize', fontsizeLAB);
set(get(gca, 'ylabel'), 'fontsize', fontsizeLAB);
set(get(gca, 'title'), 'fontsize', fontsizeLAB);
ylim([0, .8])

axpos = get(gca, 'innerpos');

box on
% export_fig(fig, ['~/Desktop/PaperSDE/Figures/OUSubs_s' num2str(10*s) '.png'], '-nocrop', '-painters', '-m5')
% % export_fig(fig, ['~/Desktop/Project/Reports/PresentationICL_20/Figures/OUSubs_s' num2str(10*s) '.png'], '-nocrop', '-painters', '-m5')


for b = 5 %[1 5]
    aFil = dlmread(['ResultsFilter_s' num2str(s,'%.1f') '_b' num2str(b) '.txt']);
    aFil_m = aFil(1:ceil(length(aSub)/2));
    aFil_s = aFil(ceil(length(aSub)/2)+1:end);

    fig = createFigure(W, H, 'enhanced', enhanced);
    
    plot([0, 1], [aHom, aHom], 'k--')
    hold on
    plot(zetas, aFil_m, 'k-o', 'markersize', 3)
%     plot(zetas, aFil_m + 2*aFil_s, 'k--')
%     plot(zetas, aFil_m - 2*aFil_s, 'k--')
    
    xlabel('$\zeta$', 'interpreter', 'laTeX')
    
    title(['Filter, $\beta = ', num2str(b), '$'], 'interpreter', 'latex')
    
    set(gca, 'fontsize', fontsizeTICK);
    set(get(gca, 'xlabel'), 'fontsize', fontsizeLAB);
    set(get(gca, 'ylabel'), 'fontsize', fontsizeLAB);
    set(get(gca, 'title'), 'fontsize', fontsizeLAB);
    ylim([0, .8])
    set(gca, 'innerpos', axpos)
    
    box on
    % export_fig(fig, ['~/Desktop/PaperSDE/Figures/OUFilt_s' num2str(10*s) '_b' num2str(b) '.png'], '-nocrop', '-painters', '-m5')
end
