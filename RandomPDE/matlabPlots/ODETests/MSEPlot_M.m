clc; clear; close all
addpath('resultsMSE_M')
%% 

% results = dlmread('outputMSE_ET.txt');
results = dlmread('outputMSE_RK4.txt');
errNum = results(1:end, 2);
M = results(1:end, 1);

%% PLOT FOR TEX

h = 0.01;
p = 4;
q = p;
err = @(M) (h.^(2*min(q, 2*p)) + 1 ./ M * h.^(2*min(q, p))).^0.5;

enhanced = 1;
fontsizeLAB = getLatexTextSize('normalsize', 'enhanced', enhanced);
fontsizeTICK = getLatexTextSize('small', 'enhanced', enhanced);

W = 6; H = 6;
fig = createFigure(W, H, 'enhanced',enhanced);

loglog(M, errNum, 'ko-', 'linewidth', .5, 'markersize', 4)
hold on
loglog(M, 30 * err(M), 'k--', 'linewidth', .5)

xlabel('$M$', 'interpreter', 'laTeX')
set(gca, 'fontsize', fontsizeTICK);
set(get(gca, 'xlabel'), 'fontsize', fontsizeLAB);
set(get(gca, 'ylabel'), 'fontsize', fontsizeLAB);
% title('ET', 'interpreter', 'latex')
title('RK4', 'interpreter', 'latex')

% print -depsc2 ../../../Reports/PaperRTSRK_18/VERSION17/Paper/MonteCarloET_M.eps
print -depsc2 ../../../Reports/PaperRTSRK_18/VERSION17/Paper/MonteCarloRK4_M.eps