clc; clear; close all
% addpath('resultsOU_epsSmall')
addpath('resultsOU_epsBig')
% addpath('resultsQU_epsBig2')

%%
s = 1;
zetas = linspace(0, 1, 11);

W = 3.5; H = 3.5;
enhanced = 3;

hom = dlmread(['ResultsHom_s' num2str(s,'%.0f') '.txt']);
aHom = hom(2);
aSub = dlmread(['ResultsSub_s' num2str(s,'%.0f') '.txt']);
aSub_m = aSub(1:ceil(length(aSub)/2));
aSub_s = aSub(ceil(length(aSub)/2)+1:end);

fig = createFigure(W, H, 'enhanced', enhanced);
fontsizeLAB = getLatexTextSize('normalsize', 'enhanced', 1);
fontsizeTICK = getLatexTextSize('small', 'enhanced', 1);

plot([0, 1], [hom(1), hom(1)], 'k-.')
hold on
plot(zetas, aSub_m, 'k')
plot(zetas, aSub_m + 2*aSub_s, 'k--')
plot(zetas, aSub_m - 2*aSub_s, 'k--')

xlabel('$\zeta$', 'interpreter', 'laTeX')

title('Subsampling', 'interpreter', 'latex')

set(gca, 'fontsize', fontsizeTICK);
set(get(gca, 'xlabel'), 'fontsize', fontsizeLAB);
set(get(gca, 'ylabel'), 'fontsize', fontsizeLAB);

yLim = get(gca, 'ylim');
set(gca, 'ylim', [0, 2]);
yLim = get(gca, 'ylim');
axpos = get(gca, 'innerpos');

box on
% export_fig(fig, ['../../../Reports/DraftMultiSDE_19/CaltechNotes_2/Figures/A2_sub_s' num2str(10*s) '.png'], '-painters', '-m5')


for b = 1:4
    
    aFil = dlmread(['ResultsFilter_s' num2str(s,'%.0f') '_b' num2str(b) '.txt']);
    aFil_m = aFil(1:ceil(length(aSub)/2));
    aFil_s = aFil(ceil(length(aSub)/2)+1:end);

    fig = createFigure(W, H, 'enhanced', enhanced);
    
    plot([0, 1], [hom(1), hom(1)], 'k-.')
    hold on
    plot(zetas, aFil_m, 'k')
    plot(zetas, aFil_m + 2*aFil_s, 'k--')
    plot(zetas, aFil_m - 2*aFil_s, 'k--')
    
    xlabel('$\zeta$', 'interpreter', 'laTeX')
    
    title(['Filter, $\beta = ', num2str(b), '$'], 'interpreter', 'latex')
    
    set(gca, 'fontsize', fontsizeTICK);
    set(get(gca, 'xlabel'), 'fontsize', fontsizeLAB);
    set(get(gca, 'ylabel'), 'fontsize', fontsizeLAB);
    set(gca, 'ylim', yLim)
    set(gca, 'innerpos', axpos)
    
    box on
%     export_fig(fig, ['../../../Reports/DraftMultiSDE_19/CaltechNotes_2/Figures/A2_filt_s' num2str(10*s) '_b' num2str(b) '.png'], '-painters', '-m5')
end
