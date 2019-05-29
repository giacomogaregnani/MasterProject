clc; clear; close all
addpath('resultsWeak')
%%

enhanced = 1;
fontsizeLAB = getLatexTextSize('normalsize', 'enhanced', enhanced);
fontsizeTICK = getLatexTextSize('small', 'enhanced', enhanced);
W = 5; H = 5;
fig = createFigure(W, H, 'enhanced',enhanced);

% q = 2;
% p = [0.5, 1, 1.5];
q = 4;
p = [1, 1.5, 2, 2.5];
markers = 'o<s+';
linestyle = {'-', '--', '-.'};

lastorder = 0;

for i = 1 : length(p)
    
%     errNum = dlmread(['outputWEAK_ET', num2str(i), '.txt']);
    errNum = dlmread(['outputWEAK_RK4', num2str(i), '.txt']);

    h = errNum(:, 1);
    err = errNum(:, 2);
    
    l(i) = loglog(errNum(:, 1), errNum(:, 2), 'ko-', 'linewidth', .5, 'marker', markers(i), 'markersize', 4);
    hold on

    order = min(2*p(i), q);
    if order > lastorder
        loglog(h, 6*h.^order, 'k', 'linestyle', linestyle{i}, 'linewidth', .5)
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
% set(gca, 'yTick', [1e-8, 1e-5, 1e-2])
set(gca, 'xTick', [1e-2, 1e-1])
set(gca, 'fontsize', fontsizeTICK);
title('ET', 'interpreter', 'latex')
% title('RK4', 'interpreter', 'latex')
set(get(gca, 'xlabel'), 'fontsize', fontsizeLAB);
set(get(gca, 'ylabel'), 'fontsize', fontsizeLAB);

% print -depsc2 ../../../Reports/PaperRTSRK_18/VERSION15/WeakET.eps
% print -depsc2 ../../../Reports/PaperRTSRK_18/VERSION15/WeakRK4.eps