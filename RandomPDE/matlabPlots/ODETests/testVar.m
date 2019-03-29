clc; clear; close all;

%%
nMC = 100;
xDet = dlmread('testVardet.txt');
x = dlmread('testVar.txt');
N = size(x, 1) / nMC;

figure
hold on
for i = 1 : nMC
    idx = [(i-1)*N + 1 : i*N];
    plot(x(idx, 1), x(idx, 2), 'color', [0.7, 0.7, 0.7]);
    xlabel('$t$', 'interpreter', 'laTeX')
end
plot(xDet(:, 1), xDet(:, 2), 'k', 'linewidth', 1)
box on


figure
stdY = [];
for i = 1 : N
    population = x(i : N : end, 2);
    stdY = [stdY, std(population)];    
end
plot(stdY)
hold on
plot(stdY(2:end) - stdY(1:end-1))
