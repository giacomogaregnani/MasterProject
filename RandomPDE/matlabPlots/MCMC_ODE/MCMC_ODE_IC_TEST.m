clc; clear; close all

resultsEE = dlmread('ODE_IC_ERREST_TEST.txt');
resultsDet = dlmread('ODE_IC_DET_TEST.txt');
resultsProb = dlmread('ODE_IC_PROB_TEST.txt');

yStar = dlmread('observationsShort.txt'); yStar = yStar(2);

figure
hold on

[f, xi] = ksdensity(resultsEE);
plot(xi, f)

[f, xi] = ksdensity(resultsDet);
plot(xi, f)

[f, xi] = ksdensity(resultsProb);
plot(xi, f)

% Truth 
lambda = -1;
t = 1;
sigmaE = 1e-3;
m =  exp(lambda * t) * yStar / (sigmaE^2 + exp(2* lambda * t));
s = sigmaE^2 / (sigmaE^2 + exp(2 * lambda * t));

x = linspace(m - 4 * sqrt(s), m + 4 * sqrt(s), 1e3);
truth = normpdf(x, m, sqrt(s));

hold on
plot(x, truth, 'k');