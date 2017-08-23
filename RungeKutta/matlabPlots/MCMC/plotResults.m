clc; clear; close all;

%% READ RESULTS

det = dlmread('samplesDet.txt');
prob = dlmread('samplesProb.txt');

% det = exp(det);
% prob = exp(prob);

nMCMC = size(det, 1);

detPlot = unique(det, 'rows');
probPlot = unique(prob, 'rows');

figure
scatter3(detPlot(:, 1), detPlot(:, 2), detPlot(:, 3), 0.5, 'r.')
hold on
scatter3(probPlot(:, 1), probPlot(:, 2), probPlot(:, 3), 0.5, 'b.')
plot3(10, 28, 8/3, 'xk', 'markersize', 10);
% plot3(0.2, 0.2, 3.0, 'xk', 'markersize', 10)
view(2)

MCMCEstDet = mean(det);
MCMCEstProb = mean(prob);


%% Create grids

xTestDet = linspace(min(det(:, 1)), max(det(:, 1)), 400);
yTestDet = linspace(min(det(:, 2)), max(det(:, 2)), 400);

xTestProb = linspace(min(prob(:, 1)), max(prob(:, 1)), 400);
yTestProb = linspace(min(prob(:, 2)), max(prob(:, 2)), 400);

[XXDet, YYDet] = meshgrid(xTestDet, yTestDet);
[XXProb, YYProb] = meshgrid(xTestProb, yTestProb);

%% Compute densities

nKDE = floor(nMCMC^(1/3));
densDet = akde(det, [det(:, 1), det(:, 2), MCMCEstDet(3)*ones(nMCMC, 1)], nKDE);
densProb = akde(prob, [prob(:, 1), prob(:, 2), MCMCEstProb(3)*ones(nMCMC, 1)], nKDE);

ZZdet = griddata(det(:, 1), det(:, 2), densDet, XXDet, YYDet);
ZZdet(isnan(ZZdet)) = 0;

ZZProb = griddata(prob(:, 1), prob(:, 2), densProb, XXProb, YYProb);
ZZProb(isnan(ZZProb)) = 0;

%% Plot results

figure
contour(XXDet, YYDet, ZZdet, 12, 'b');
% xLim = get(gca, 'xlim');
% yLim = get(gca, 'ylim');
hold on
plot(0.2, 0.2, 'xk', 'LineWidth', 3)
axis tight
hold on
contour(XXProb, YYProb, ZZProb, 12, 'r');

% scatter(prob(:, 1), prob(:, 2), 0.5, 'r.')
