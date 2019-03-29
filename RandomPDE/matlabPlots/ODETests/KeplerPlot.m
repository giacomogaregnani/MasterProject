clc; clear; close all

%% Read data
xTotAdd = dlmread('IMFullAdd.txt');
xTotStep = dlmread('IMFullStep.txt');
h = xTotAdd(2, 1) - xTotAdd(1, 1);

%% Figure param
W = 4; H = 4;
enhanced = 1;
darkgrey = [0.3, 0.3, 0.3];
lightgrey = [0.6, 0.6, 0.6];
offset = 0.1;

%% Step noise

fig = createFigure(W, H, 'enhanced',enhanced);
fontsizeLAB = getLatexTextSize('normalsize', 'enhanced', enhanced);
fontsizeTICK = getLatexTextSize('small', 'enhanced', enhanced);

plot(xTotStep(1:200/h, 4), xTotStep(1:200/h, 5), 'color', darkgrey);

set(gca, 'fontsize', fontsizeTICK);
set(get(gca, 'xlabel'), 'fontsize', fontsizeLAB);
set(get(gca, 'ylabel'), 'fontsize', fontsizeLAB);
set(gca, 'ytick', [-1 1])
% set(gca, 'yticklabel', [])
set(gca, 'xtick', [-1 1])
% set(gca, 'xticklabel', [])

xlabel('$w_1$', 'interpreter', 'latex')
ylabel('$w_2$', 'interpreter', 'latex')

axis equal
ax = axis;
% text(ax(1) + offset, -1, ['-1'],...
%     'HorizontalAlignment','Left', 'fontsize', fontsizeTICK);
% text(ax(1) + offset, 1, ['1'],...
%     'HorizontalAlignment','Left', 'fontname', 'CMUSerif', 'fontsize', fontsizeTICK);
% text(-1, ax(3) + offset, -1, ['-1'],...
%     'HorizontalAlignment','Center', 'fontname', 'CMUSerif', 'fontsize', fontsizeTICK);
% text(1, ax(3) + offset, 1, ['1'],...
%     'HorizontalAlignment','Center', 'fontname', 'CMUSerif', 'fontsize', fontsizeTICK);
xLim = get(gca, 'xlim');
yLim = get(gca, 'ylim');


print -depsc2 ../../../Reports/PaperRTSRK_18/VERSION12/KeplerOne.eps

fig = createFigure(W, H, 'enhanced',enhanced);
fontsizeLAB = getLatexTextSize('normalsize', 'enhanced', enhanced);
fontsizeTICK = getLatexTextSize('small', 'enhanced', enhanced);

plot(xTotStep(3800/h:4000/h, 4), xTotStep(3800/h:4000/h, 5), 'color', darkgrey);

set(gca, 'fontsize', fontsizeTICK);
set(get(gca, 'xlabel'), 'fontsize', fontsizeLAB);
set(get(gca, 'ylabel'), 'fontsize', fontsizeLAB);
set(gca, 'ytick', [-1 1])
% set(gca, 'yticklabel', [])
set(gca, 'xtick', [-1 1])
% set(gca, 'xticklabel', [])

xlabel('$w_1$', 'interpreter', 'latex')
ylabel('$w_2$', 'interpreter', 'latex')

axis equal
ax = axis;
% text(ax(1) + offset, -1, ['-1'],...
%     'HorizontalAlignment','Left', 'fontname', 'CMUSerif', 'fontsize', fontsizeTICK);
% text(ax(1) + offset, 1, ['1'],...
%     'HorizontalAlignment','Left', 'fontname', 'CMUSerif', 'fontsize', fontsizeTICK);
% text(-1, ax(3) + offset, -1, ['-1'],...
%     'HorizontalAlignment','Center', 'fontname', 'CMUSerif', 'fontsize', fontsizeTICK);
% text(1, ax(3) + offset, 1, ['1'],...
%     'HorizontalAlignment','Center', 'fontname', 'CMUSerif', 'fontsize', fontsizeTICK);

print -depsc2 ../../../Reports/PaperRTSRK_18/VERSION12/KeplerTwo.eps

%% Additive noise

fig = createFigure(W, H, 'enhanced',enhanced);
fontsizeLAB = getLatexTextSize('normalsize', 'enhanced', enhanced);
fontsizeTICK = getLatexTextSize('small', 'enhanced', enhanced);

plot(xTotAdd(1:200/h, 4), xTotAdd(1:200/h, 5), 'color', darkgrey);

set(gca, 'fontsize', fontsizeTICK);
set(get(gca, 'xlabel'), 'fontsize', fontsizeLAB);
set(get(gca, 'ylabel'), 'fontsize', fontsizeLAB);
set(gca, 'ytick', [-1 1])
% set(gca, 'yticklabel', [])
set(gca, 'xtick', [-1 1])
% set(gca, 'xticklabel', [])

axis equal
set(gca, 'xlim', xLim)
set(gca, 'ylim', yLim)
ax = axis;
% text(ax(1) + offset, -1, ['-1'],...
%     'HorizontalAlignment','Left', 'fontname', 'CMUSerif', 'fontsize', fontsizeTICK);
% text(ax(1) + offset, 1, ['1'],...
%     'HorizontalAlignment','Left', 'fontname', 'CMUSerif', 'fontsize', fontsizeTICK);
% text(-1, ax(3) + offset, -1, ['-1'],...
%     'HorizontalAlignment','Center', 'fontname', 'CMUSerif', 'fontsize', fontsizeTICK);
% text(1, ax(3) + offset, 1, ['1'],...
%     'HorizontalAlignment','Center', 'fontname', 'CMUSerif', 'fontsize', fontsizeTICK);

xlabel('$w_1$', 'interpreter', 'latex')
ylabel('$w_2$', 'interpreter', 'latex')

print -depsc2 ../../../Reports/PaperRTSRK_18/VERSION12/KeplerOneAdd.eps

fig = createFigure(W, H, 'enhanced',enhanced);
fontsizeLAB = getLatexTextSize('normalsize', 'enhanced', enhanced);
fontsizeTICK = getLatexTextSize('small', 'enhanced', enhanced);

plot(xTotAdd(200/h:400/h, 4), xTotAdd(200/h:400/h, 5), 'color', darkgrey);

set(gca, 'fontsize', fontsizeTICK);
set(get(gca, 'xlabel'), 'fontsize', fontsizeLAB);
set(get(gca, 'ylabel'), 'fontsize', fontsizeLAB);
set(gca, 'ytick', [-1 1])
% set(gca, 'yticklabel', [])
set(gca, 'xtick', [-1 1])
% set(gca, 'xticklabel', [])

axis equal
set(gca, 'xlim', xLim)
set(gca, 'ylim', yLim)
ax = axis;
% text(ax(1) + offset, -1, ['-1'],...
%     'HorizontalAlignment','Left', 'fontname', 'CMUSerif', 'fontsize', fontsizeTICK);
% text(ax(1) + offset, 1, ['1'],...
%     'HorizontalAlignment','Left', 'fontname', 'CMUSerif', 'fontsize', fontsizeTICK);
% text(-1, ax(3) + offset, -1, ['-1'],...
%     'HorizontalAlignment','Center', 'fontname', 'CMUSerif', 'fontsize', fontsizeTICK);
% text(1, ax(3) + offset, 1, ['1'],...
%     'HorizontalAlignment','Center', 'fontname', 'CMUSerif', 'fontsize', fontsizeTICK);

xlabel('$w_1$', 'interpreter', 'latex')
ylabel('$w_2$', 'interpreter', 'latex')

print -depsc2 ../../../Reports/PaperRTSRK_18/VERSION12/KeplerTwoAdd.eps

%% Momentum

H = 4;
W = 12;

fig = createFigure(W, H, 'enhanced',enhanced);
fontsizeLAB = getLatexTextSize('normalsize', 'enhanced', enhanced);
fontsizeTICK = getLatexTextSize('small', 'enhanced', enhanced);

t = xTotAdd(2:end, 1);
semilogy(t, abs(xTotAdd(2:end, 6) - xTotAdd(1, 6)), 'color', lightgrey);
hold on
semilogy(t, abs(xTotStep(2:end, 6) - xTotStep(1, 6)), 'color', 'k');

xlabel('$t$', 'interpreter', 'LaTex');
ylabel('$|I(v, w) - I(v_0, w_0)|$', 'interpreter', 'LaTex');
set(gca, 'fontsize', fontsizeTICK);
set(get(gca, 'xlabel'), 'fontsize', fontsizeLAB);
set(get(gca, 'ylabel'), 'fontsize', fontsizeLAB);
xlim([0, t(end)])
set(gca, 'xtick', 0:500:4000);

legend({'Additive noise', 'Random time step'}, 'interpreter', 'latex', 'Location', 'E')

print -depsc2 ../../../Reports/PaperRTSRK_18/VERSION12/KeplerMom.eps




