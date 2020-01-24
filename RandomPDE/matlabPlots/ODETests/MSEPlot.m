clc; clear; close all
addpath('resultsMSE')
%% 
% For the parameters see the paper with VERSION >= 15

h = 0.125 ./ 2.^[0:7];

errNum = dlmread(['outputMSE_ET.txt']);
errNum = errNum(2:end, :);

% M = 1e3;
% p = 1; q = 2;
M = 1e4;
p = 1.5; q = 4;

err = @(h) (h.^(2*min(q, 2*p)) + 1 / M * h.^(2*min(q, p))).^0.5;

% loglog(h, err(h), 'k--', 'linewidth', .2)

%% PLOT FOR TEX

enhanced = 1;
fontsizeLAB = getLatexTextSize('normalsize', 'enhanced', enhanced);
fontsizeTICK = getLatexTextSize('small', 'enhanced', enhanced);

W = 6; H = 6;
fig = createFigure(W, H, 'enhanced',enhanced);

loglog(errNum(:, 1), errNum(:, 2), 'ko-', 'linewidth', .5, 'markersize', 4)
hold on
loglog(h, err(h), 'k--', 'linewidth', .5)

% legend('result', 'estimate', 'location', 'best')
xlabel('$h$', 'interpreter', 'laTeX')
% ylabel('$\mathrm{MSE}^{1/2}$', 'interpreter', 'laTeX')
xlim([1e-4, 1e0])
% ylim([1e-6, 1e0])
set(gca, 'yTick', [1e-6, 1e-3, 1e0]) 
% set(gca, 'yTick', [1e-4, 1e-2, 1e0]) 
set(gca, 'fontsize', fontsizeTICK);
set(get(gca, 'xlabel'), 'fontsize', fontsizeLAB);
set(get(gca, 'ylabel'), 'fontsize', fontsizeLAB);
% title('ET', 'interpreter', 'latex')
title('RK4', 'interpreter', 'latex')

% print -depsc2 ../../../Reports/PaperRTSRK_18/VERSION15/MonteCarloET.eps
% print -depsc2 ../../../Reports/PaperRTSRK_18/VERSION15/MonteCarloRK4.eps