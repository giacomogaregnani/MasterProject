clc; clear; close all

u = dlmread('MCMC_PDE1.txt');

meanU = mean(u);
stdU = std(u);

x = linspace(0, 1, 1001);
xEx = linspace(0, 1, 1000);

plot(x, meanU, 'b')
hold on
plot(x, meanU + 1.96 * stdU, 'b--')
plot(x, meanU - 1.96 * stdU, 'b--')
plot(xEx, 1 + 0.5 * sin(2 * pi * xEx), 'k')
