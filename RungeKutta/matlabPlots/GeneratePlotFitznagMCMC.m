clear; clc; close all;
%% Figure parameters

W = 6; H = 6;
enhanced = 1;
fontsizeLAB = getLatexTextSize('normalsize', 'enhanced', enhanced);
fontsizeTICK = getLatexTextSize('small', 'enhanced', enhanced);

%% Load results "first row"

load('resultsMCMCh01')

%% Generate plots "first row"

fig = createFigure(W, H, 'enhanced',enhanced);

contourf(XX, YY, ZZStep, 12);
hold on
plot(0.2, 0.2, 'r.', 'MarkerSize', 20)

set(gca, 'fontsize', fontsizeTICK);

axis square

print('../../Reports/RandTimeStep/VERSION4/MCMCFitznagStep01.png', '-r2000', '-dpng')

fig = createFigure(W, H, 'enhanced',enhanced);

contourf(XX, YY, ZZAdd, 12);
hold on
plot(0.2, 0.2, 'r.', 'MarkerSize', 20)

set(gca, 'fontsize', fontsizeTICK);

axis square

print('../../Reports/RandTimeStep/VERSION4/MCMCFitznagAdd01.png', '-r2000', '-dpng')


%% Load results "second row"

load('resultsMCMCh001')

%% Generate plots "second row"

fig = createFigure(W, H, 'enhanced',enhanced);

contourf(XX, YY, ZZStep, 12);
hold on
plot(0.2, 0.2, 'r.', 'MarkerSize', 20)

set(gca, 'fontsize', fontsizeTICK);

axis square

print('../../Reports/RandTimeStep/VERSION4/MCMCFitznagStep001.png', '-r2000', '-dpng')

fig = createFigure(W, H, 'enhanced',enhanced);

contourf(XX, YY, ZZAdd, 12);
hold on
plot(0.2, 0.2, 'r.', 'MarkerSize', 20)

set(gca, 'fontsize', fontsizeTICK);

axis square

print('../../Reports/RandTimeStep/VERSION4/MCMCFitznagAdd001.png', '-r2000', '-dpng')


