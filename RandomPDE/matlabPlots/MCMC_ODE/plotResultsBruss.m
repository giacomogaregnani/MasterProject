clc; clear; close all
%%

results = dlmread('MCMC_BRUSS1.txt');
thetas = dlmread('MCMC_BRUSS1theta.txt');

resultsProb = dlmread('MCMC_BRUSS_PROB1.txt');
thetasProb =  dlmread('MCMC_BRUSS_PROB1theta.txt');

%%

meanMCMC = mean(results);
varMCMC = var(results);

meanMCMCProb = mean(resultsProb);
varMCMCProb = var(resultsProb);

X = linspace(0, 1, 32);
%%

figure
hold on

plot(X, meanMCMC, 'b')
plot(X, meanMCMC + 2*sqrt(varMCMC), 'b--')
plot(X, meanMCMC - 2*sqrt(varMCMC), 'b--')

plot(X, meanMCMCProb, 'r')
plot(X, meanMCMCProb + 2*sqrt(varMCMCProb), 'r--')
plot(X, meanMCMCProb - 2*sqrt(varMCMCProb), 'r--')

xEx = linspace(0, 1, 1000);
plot(xEx, 1 + sqrt(2) / (2 * pi) * sin(2 * pi * xEx), 'k')

figure
hold on

[f, xi] = ksdensity(thetas(:, 2));
plot(xi, f)

[fProb, xiProb] = ksdensity(thetasProb(:, 2));
plot(xiProb, fProb)