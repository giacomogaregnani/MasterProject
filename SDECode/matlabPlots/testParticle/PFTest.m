% clc; clear; close all
%%

obs = dlmread('testSol.txt');
PF = dlmread('test.txt');

t = linspace(0, 10, length(obs));

figure
plot(t, PF(1:end-1,:), 'color', 0.6 * ones(1, 3))
hold on
plot(t, PF(end,:), 'r')
plot(t, obs, 'k')

figure
l = dlmread('testLik.txt');
hold on
H = histogram(l, 'normalization', 'pdf');
hold on
meanl = mean(l);
plot([meanl, meanl], [0, max(H.Values)], 'k--', 'linewidth', 2)
s = std(l)