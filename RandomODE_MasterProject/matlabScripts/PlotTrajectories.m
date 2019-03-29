clc; clear; close all
%% READ DATA

data = dlmread('output.txt');

nTrajectories = 20;

%% Plot
h = 0.05;
T = 20;
N = T / h + 1;
timeVec = 0 : h : T;

figure
ax1 = axes;
hold on
title('first component')
figure
ax2 = axes;
hold on
title('second component')

idx = N + 1;
for i = 1 : nTrajectories
    yProb = data(idx:idx+N-1, :);
    plot(ax1, timeVec, yProb(:, 1), 'color', [0.7, 0.7, 0.7])
    plot(ax2, timeVec, yProb(:, 2), 'color', [0.7, 0.7, 0.7])
    idx = idx + N;
end

yDet = data(1:N, :);
plot(ax1, timeVec, yDet(:, 1), 'k', 'LineWidth', 2)
plot(ax2, timeVec, yDet(:, 2), 'k', 'LineWidth', 2)
