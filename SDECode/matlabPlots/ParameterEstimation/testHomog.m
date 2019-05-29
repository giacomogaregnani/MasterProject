clc; clear; close all
%%

results = dlmread('TestHomog.txt');

x = results(:, 1);
xH = results(:, 2);

figure
[f, xi] = ksdensity(x);
plot(xi, f);
hold on
[f, xi] = ksdensity(xH);
plot(xi, f);
% histogram(x, 'normalization', 'pdf')
% hold on
% histogram(xH, 'normalization', 'pdf')