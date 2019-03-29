clc; clear; close all
addpath('resultsMSE')
%%

h = 0.25 ./ 2.^[1:8];

M = 1e3;
errNum = dlmread(['outputMSE1e3_ET.txt']);
errNum = errNum(2:end, :);

p = 1.5; q = 2;

err = @(h) (h.^(2*min(q, 2*p-1)) + 1 / M * h.^(2*min(q, p-0.5))).^0.5;

err1 = @(h) h.^(min(q, 2*p-1));
err2 = @(h) 1 / sqrt(M) * h.^(min(q, p-0.5));

% loglog(h, err(h), 'k--', 'linewidth', .2)

%% PLOT FOR TEX

enhanced = 1;
fontsizeLAB = getLatexTextSize('normalsize', 'enhanced', enhanced);
fontsizeTICK = getLatexTextSize('small', 'enhanced', enhanced);

W = 4.2; H = 4;
fig = createFigure(W, H, 'enhanced',enhanced);

loglog(errNum(:, 1), errNum(:, 2), 'ko-', 'linewidth', .5, 'markersize', 4)
hold on
loglog(h, err(h), 'k--', 'linewidth', .5)

% legend('result', 'estimate', 'location', 'best')
xlabel('$h$', 'interpreter', 'laTeX')
ylabel('$\mathrm{MSE}^{1/2}$', 'interpreter', 'laTeX')
xlim([1e-4, 1e0])
ylim([1e-6, 1e0])
set(gca, 'yTick', [1e-6, 1e-3, 1e0]) 
set(gca, 'fontsize', fontsizeTICK);
set(get(gca, 'xlabel'), 'fontsize', fontsizeLAB);
set(get(gca, 'ylabel'), 'fontsize', fontsizeLAB);

print -depsc2 ../../../Reports/PaperRTSRK_18/VERSION11/MonteCarloET.eps