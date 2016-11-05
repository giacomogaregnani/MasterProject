clear
close all

results = dlmread('../data/resultsStab.txt');

[twoN, size] = size(results);
N = twoN / 2;

resultsStab = results(1 : N, :);
resultsNStab = results(N+1 : end, :);

t = 0 : 0.1 : 20;

plot(t, resultsStab(:, 1))
hold on
plot(t, resultsNStab(:, 1), 'r')