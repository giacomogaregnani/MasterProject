clc; clear; close all

yProb = dlmread('FitzNag.txt');
nMC = 50;
nTimes = size(yProb, 1) / nMC;

W = 10; H = 3; 
fontsizeLAB = getLatexTextSize('normalsize', 'enhanced', 1);
fontsizeTICK = getLatexTextSize('footnotesize', 'enhanced', 1);

createFigure(W, H, 'enhanced', 1);
hold on
for i = 1 : 20
    y = yProb((i-1)*nTimes+1:i*nTimes, :);
    plot(y(:, 1), y(:, 2), 'color', 0.7 * ones(3, 1));
end

yDet = dlmread('FitzNagdet.txt');
plot(yDet(:, 1), yDet(:, 2), 'k')
xlabel('$t$', 'interpreter', 'laTeX', 'fontsize', fontsizeTICK)
ylabel('$y_1$', 'interpreter', 'laTeX', 'fontsize', fontsizeTICK)
ylim([-4, 4])


box on
% print -depsc2 ../../../Reports/PresentationMaxPlanck/FitzNag.eps

createFigure(W, H, 'enhanced', 1);
hold on
for i = 1 : 20
    y = yProb((i-1)*nTimes+1:i*nTimes, :);
    plot(y(:, 1), y(:, 3), 'color', 0.7 * ones(3, 1));
end

yDet = dlmread('FitzNagdet.txt');
plot(yDet(:, 1), yDet(:, 3), 'k')
xlabel('$t$', 'interpreter', 'laTeX', 'fontsize', fontsizeTICK)
ylabel('$y_2$', 'interpreter', 'laTeX', 'fontsize', fontsizeTICK)
ylim([-2, 2])

box on
% print -depsc2 ../../../Reports/PresentationMaxPlanck/FitzNag2.eps
