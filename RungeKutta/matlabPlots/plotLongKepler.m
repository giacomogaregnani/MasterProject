clc; clear; close all

%% Read data
xTotAdd = dlmread('IMFullAdd.txt');
xTotStep = dlmread('IMFullStep.txt');

%% Figure param
W = 6; H = 6;
enhanced = 1;
darkgrey = [0.3, 0.3, 0.3];
lightgrey = [0.6, 0.6, 0.6];
offset = 0.1;


%% Step noise

fig = createFigure(W, H, 'enhanced',enhanced);
fontsizeLAB = getLatexTextSize('normalsize', 'enhanced', enhanced);
fontsizeTICK = getLatexTextSize('small', 'enhanced', enhanced);

plot(xTotStep(1:20000, 3), xTotStep(1:20000, 4), 'color', darkgrey);

set(gca, 'fontsize', fontsizeTICK);
set(get(gca, 'xlabel'), 'fontsize', fontsizeLAB);
set(get(gca, 'ylabel'), 'fontsize', fontsizeLAB);
set(gca, 'ytick', [-1 1])
set(gca, 'yticklabel', [])
set(gca, 'xtick', [-1 1])
set(gca, 'xticklabel', [])
axis equal
ax = axis;
text(ax(1) + offset, -1, ['-1'],...
    'HorizontalAlignment','Left', 'fontsize', fontsizeTICK);
text(ax(1) + offset, 1, ['1'],...
    'HorizontalAlignment','Left', 'fontname', 'CMUSerif', 'fontsize', fontsizeTICK);
text(-1, ax(3) + offset, -1, ['-1'],...
    'HorizontalAlignment','Center', 'fontname', 'CMUSerif', 'fontsize', fontsizeTICK);
text(1, ax(3) + offset, 1, ['1'],...
    'HorizontalAlignment','Center', 'fontname', 'CMUSerif', 'fontsize', fontsizeTICK);
xLim = get(gca, 'xlim');
yLim = get(gca, 'ylim');

print -depsc2 ../../Reports/RandTimeStep/VERSION3/KeplerOne.eps

fig = createFigure(W, H, 'enhanced',enhanced);
fontsizeLAB = getLatexTextSize('normalsize', 'enhanced', enhanced);
fontsizeTICK = getLatexTextSize('small', 'enhanced', enhanced);

plot(xTotStep(380000:400000, 3), xTotStep(380000:400000, 4), 'color', darkgrey);

set(gca, 'fontsize', fontsizeTICK);
set(get(gca, 'xlabel'), 'fontsize', fontsizeLAB);
set(get(gca, 'ylabel'), 'fontsize', fontsizeLAB);
set(gca, 'ytick', [-1 1])
set(gca, 'yticklabel', [])
set(gca, 'xtick', [-1 1])
set(gca, 'xticklabel', [])

axis equal
ax = axis;
text(ax(1) + offset, -1, ['-1'],...
    'HorizontalAlignment','Left', 'fontname', 'CMUSerif', 'fontsize', fontsizeTICK);
text(ax(1) + offset, 1, ['1'],...
    'HorizontalAlignment','Left', 'fontname', 'CMUSerif', 'fontsize', fontsizeTICK);
text(-1, ax(3) + offset, -1, ['-1'],...
    'HorizontalAlignment','Center', 'fontname', 'CMUSerif', 'fontsize', fontsizeTICK);
text(1, ax(3) + offset, 1, ['1'],...
    'HorizontalAlignment','Center', 'fontname', 'CMUSerif', 'fontsize', fontsizeTICK);

print -depsc2 ../../Reports/RandTimeStep/VERSION3/KeplerTwo.eps

%% Additive noise

fig = createFigure(W, H, 'enhanced',enhanced);
fontsizeLAB = getLatexTextSize('normalsize', 'enhanced', enhanced);
fontsizeTICK = getLatexTextSize('small', 'enhanced', enhanced);

plot(xTotAdd(1:20000, 3), xTotAdd(1:20000, 4), 'color', darkgrey);

set(gca, 'fontsize', fontsizeTICK);
set(get(gca, 'xlabel'), 'fontsize', fontsizeLAB);
set(get(gca, 'ylabel'), 'fontsize', fontsizeLAB);
set(gca, 'ytick', [-1 1])
set(gca, 'yticklabel', [])
set(gca, 'xtick', [-1 1])
set(gca, 'xticklabel', [])

axis equal
set(gca, 'xlim', xLim)
set(gca, 'ylim', yLim)
ax = axis;
text(ax(1) + offset, -1, ['-1'],...
    'HorizontalAlignment','Left', 'fontname', 'CMUSerif', 'fontsize', fontsizeTICK);
text(ax(1) + offset, 1, ['1'],...
    'HorizontalAlignment','Left', 'fontname', 'CMUSerif', 'fontsize', fontsizeTICK);
text(-1, ax(3) + offset, -1, ['-1'],...
    'HorizontalAlignment','Center', 'fontname', 'CMUSerif', 'fontsize', fontsizeTICK);
text(1, ax(3) + offset, 1, ['1'],...
    'HorizontalAlignment','Center', 'fontname', 'CMUSerif', 'fontsize', fontsizeTICK);

print -depsc2 ../../Reports/RandTimeStep/VERSION3/KeplerOneAdd.eps

fig = createFigure(W, H, 'enhanced',enhanced);
fontsizeLAB = getLatexTextSize('normalsize', 'enhanced', enhanced);
fontsizeTICK = getLatexTextSize('small', 'enhanced', enhanced);

plot(xTotAdd(20000:40000, 3), xTotAdd(20000:40000, 4), 'color', darkgrey);

set(gca, 'fontsize', fontsizeTICK);
set(get(gca, 'xlabel'), 'fontsize', fontsizeLAB);
set(get(gca, 'ylabel'), 'fontsize', fontsizeLAB);
set(gca, 'ytick', [-1 1])
set(gca, 'yticklabel', [])
set(gca, 'xtick', [-1 1])
set(gca, 'xticklabel', [])

axis equal
set(gca, 'xlim', xLim)
set(gca, 'ylim', yLim)
ax = axis;
text(ax(1) + offset, -1, ['-1'],...
    'HorizontalAlignment','Left', 'fontname', 'CMUSerif', 'fontsize', fontsizeTICK);
text(ax(1) + offset, 1, ['1'],...
    'HorizontalAlignment','Left', 'fontname', 'CMUSerif', 'fontsize', fontsizeTICK);
text(-1, ax(3) + offset, -1, ['-1'],...
    'HorizontalAlignment','Center', 'fontname', 'CMUSerif', 'fontsize', fontsizeTICK);
text(1, ax(3) + offset, 1, ['1'],...
    'HorizontalAlignment','Center', 'fontname', 'CMUSerif', 'fontsize', fontsizeTICK);

print -depsc2 ../../Reports/RandTimeStep/VERSION3/KeplerTwoAdd.eps

%% Momentum

W = 14;

fig = createFigure(W, H, 'enhanced',enhanced);
fontsizeLAB = getLatexTextSize('normalsize', 'enhanced', enhanced);
fontsizeTICK = getLatexTextSize('small', 'enhanced', enhanced);

t = linspace(0.01, 0.01 * 408408, 408407);
semilogy(t, abs(xTotAdd(2:end, 5) - xTotAdd(1, 5)), 'color', lightgrey);
hold on
semilogy(t, abs(xTotStep(2:end, 5) - xTotStep(1, 5)), 'color', 'k');

xlabel('$t$', 'interpreter', 'LaTex');
ylabel('$|I(p, q) - I(p_0, q_0)|$', 'interpreter', 'LaTex');
set(gca, 'fontsize', fontsizeTICK);
set(get(gca, 'xlabel'), 'fontsize', fontsizeLAB);
set(get(gca, 'ylabel'), 'fontsize', fontsizeLAB);
xlim([0, t(end)])

legend('Additive noise', 'Random time step', 'Location', 'best')

print -depsc2 ../../Reports/RandTimeStep/VERSION3/KeplerMom.eps




