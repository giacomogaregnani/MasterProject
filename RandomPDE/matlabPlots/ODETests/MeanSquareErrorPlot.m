clc; clear; close all
addpath('resultsMeanSquare')
%%

enhanced = 1;
fontsizeLAB = getLatexTextSize('normalsize', 'enhanced', enhanced);
fontsizeTICK = getLatexTextSize('small', 'enhanced', enhanced);
W = 6; H = 6;
fig = createFigure(W, H, 'enhanced',enhanced);

p = [1, 2, 3];
q = 2;
% q = 4;
% p = [2, 3, 4, 5];
markers = 'o<s+';
linestyle = {'-', '--', '-.'};

lastorder = 0;

for i = 1 : length(p)
    
    errNum = dlmread(['outputMSC_ET', num2str(i), '.txt']);
%     errNum = dlmread(['outputMSC_RK4', num2str(i), '.txt']);

    h = errNum(:, 1);
    err = errNum(:, 2);
    
    l(i) = loglog(errNum(:, 1), errNum(:, 2), 'ko-', 'linewidth', .5, 'marker', markers(i), 'markersize', 4);
    hold on

    order = min(p(i), q);
    if order > lastorder
        loglog(h, h.^order, 'k', 'linestyle', linestyle{i}, 'linewidth', .5)
    end
    lastorder = order;
end


%% PLOT FOR TEX

xlabel('$h$', 'interpreter', 'laTeX')
% legend([l(1), l(2), l(3), l(4)], {['$p = ', num2str(p(1)), '$'], ['$p = ', num2str(p(2)), '$'], ...
%                                   ['$p = ', num2str(p(3)), '$'], ['$p = ', num2str(p(4)), '$']}, 'interpreter', 'latex', 'location', 'se')
legend([l(1), l(2), l(3)], {['$p = ', num2str(p(1)), '$'], ['$p = ', num2str(p(2)), '$'], ...
                            ['$p = ', num2str(p(3)), '$']}, 'interpreter', 'latex', 'location', 'se')
set(gca, 'yTick', [1e-4, 1e-2, 1e0])
% set(gca, 'yTick', [1e-9, 1e-6, 1e-3])
set(gca, 'xTick', [1e-2, 1e-1])
set(gca, 'fontsize', fontsizeTICK);
title('ET', 'interpreter', 'latex')
% title('RK4', 'interpreter', 'latex')
set(get(gca, 'xlabel'), 'fontsize', fontsizeLAB);
set(get(gca, 'ylabel'), 'fontsize', fontsizeLAB);

print -depsc2 ../../../Reports/PaperRTSRK_18/VERSION15/MeanSquareET.eps
% print -depsc2 ../../../Reports/PaperRTSRK_18/VERSION15/MeanSquareRK4.eps