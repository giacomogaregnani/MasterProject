clc; clear; close all
%% Read data
xStep = dlmread('resultsMCMCStepLorenz.txt');
xAdd = dlmread('resultsMCMCAdd01.txt');

%% Exact posterior

densStep = akde(xStep, [xStep(:, 1), xStep(:, 2), 8 / 3 * ones(50001, 1)]);
densAdd = akde(xAdd, [xAdd(:, 1), xAdd(:, 2), 3 * ones(50001, 1)]);


%% Create grids

xTest = linspace(0.15, 0.25, 400);
yTest = linspace(-0.1, 0.4, 400);

[XX, YY] = meshgrid(xTest, yTest);

ZZAdd = griddata(xAdd(:, 1), xAdd(:, 2), densAdd, XX, YY);
ZZAdd(isnan(ZZAdd)) = 0;

ZZStep = griddata(xStep(:, 1), xStep(:, 2), densStep, XX, YY);
ZZStep(isnan(ZZStep)) = 0;


%% Plot results

contourf(XX, YY, ZZAdd, 12);
xLim = get(gca, 'xlim');
yLim = get(gca, 'ylim');
hold on
plot(0.2, 0.2, 'xr', 'LineWidth', 3)
axis tight

figure
contourf(XX, YY, ZZStep, 12);
set(gca, 'xlim', xLim)
set(gca, 'ylim', yLim)
hold on
plot(0.2, 0.2, 'xr', 'LineWidth', 3)




