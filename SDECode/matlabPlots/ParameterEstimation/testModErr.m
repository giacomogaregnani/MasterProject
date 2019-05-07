clc; clear; close all
%% 

results = dlmread('testModErr.txt');

T = 1;
means = results(1, :);
stddevs = results(2, :);
meansPF = results(3, :);
stddevsPF = results(4, :);

t = linspace(0, 1, length(means));
plot(t, means, 'k')
hold on
plot(t, means + 2 * stddevs, 'k--')
plot(t, means - 2 * stddevs, 'k--')

plot(t, meansPF, 'r')
plot(t, meansPF + 2 * stddevs, 'r--')
plot(t, meansPF - 2 * stddevs, 'r--')

