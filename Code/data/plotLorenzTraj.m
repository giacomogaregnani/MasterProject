clc; close all;

results = dlmread('LorenzTrajSet.txt');

nSteps = size(results, 1) / 11;
t = [0; results(1 : nSteps, 1)];

ax = zeros(3,1);
fig = zeros(3,1);

y0 = [-10, -1, 40];

for k = 1 : 3
    fig(k) = figure;
    hold all
    ax(k) = gca;
end

colour = [0.65 0.65 0.65];
count = nSteps + 1;
for i = 2 : 10
    y = [y0; results(count:count+nSteps-1, 2 : end)];
    for k = 1 : 3
        plot(ax(k), t, y(:, k), 'g', 'LineWidth', 0.01, 'color', colour);
    end
    count = count + nSteps;
end

for k = 1 : 3
    yRef =[y0(k); results(1 : nSteps, k + 1)];
    plot(ax(k), t, yRef, 'k', 'LineWidth', 2);
end
