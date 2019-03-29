clc; clear; close all

yIE = dlmread('KeplerIE.txt');

fontsizeLAB = getLatexTextSize('normalsize', 'enhanced', 1);
fontsizeTICK = getLatexTextSize('footnotesize', 'enhanced', 1);
W = 5; H = 5;

createFigure(W, H, 'enhanced', 1);
plot(yIE(:, 4), yIE(:, 5), 'color', 0.7 * ones(3, 1))
set(gca, 'xLim', [-2, 2]);
set(gca, 'yLim', [-1, 2.5]);
xLim = get(gca, 'xLim');
yLim = get(gca, 'yLim');
text(-1.5, -0.75, '$h = 0.005$', 'interpreter', 'latex', 'fontsize', fontsizeTICK)
title('Symplectic Euler', 'interpreter', 'latex', 'fontsize', fontsizeTICK)
xlabel('$q_1$', 'interpreter', 'laTeX', 'fontsize', fontsizeTICK)
ylabel('$q_2$', 'interpreter', 'laTeX', 'fontsize', fontsizeTICK)

% print -depsc2 ../../../Reports/PresentationMaxPlanck/KeplerSE.eps

yEE = dlmread('KeplerEE.txt');

createFigure(W, H, 'enhanced', 1);
plot(yEE(:, 4), yEE(:, 5), 'color', 0.7 * ones(3, 1))
set(gca, 'xLim', xLim)
set(gca, 'yLim', yLim)
text(-1.5, -0.75, '$h = 0.0005$', 'interpreter', 'latex', 'fontsize', fontsizeTICK)
title('Explicit Euler', 'interpreter', 'latex', 'fontsize', fontsizeTICK)
xlabel('$q_1$', 'interpreter', 'laTeX', 'fontsize', fontsizeTICK)
ylabel('$q_2$', 'interpreter', 'laTeX', 'fontsize', fontsizeTICK)

% print -depsc2 ../../../Reports/PresentationMaxPlanck/KeplerEE.eps

%% PLOT HAMILTONIANS

Q = @(p, q) 0.5 * (p(1)^2 + p(2)^2) - 1 / norm(q);
Q0 = Q(yIE(1, 2:3), yIE(1, 4:5));

createFigure(2*W + 1, H, 'enhanced', 1);
for i = 1 : size(yIE, 1)
    y = yIE(i, :);
    HamIE(i) = Q(y(2:3), y(4:5)) - Q0;
end
plot(yIE(:, 1), HamIE, 'k')

for i = 1 : size(yEE, 1)
    y = yEE(i, :);
    HamEE(i) = Q(y(2:3), y(4:5)) - Q0;
end
hold on
plot(yEE(:, 1), HamEE, 'color', 0.7 * ones(3, 1))

title('Error on energy', 'interpreter', 'latex', 'fontsize', fontsizeTICK)
xlabel('$t$', 'interpreter', 'laTeX', 'fontsize', fontsizeTICK)
lgd = legend({'S.E., $h = 0.005$', 'E.E., $h = 0.0005$'}, 'location', 'NW', 'interpreter', 'latex');
lgd.FontSize = fontsizeTICK;

% print -depsc2 ../../../Reports/PresentationMaxPlanck/KeplerEnergy.eps





