clc; clear; close all;
%%

x = dlmread('MCMC_ODE1.txt');
thetaEx = [0.2, 0.2, 3.0];

for i = 1 : length(thetaEx)
    figure
    [f, xi] = ksdensity(exp(x(:, i)));
    plot(xi, f)    
    hold on
    plot([thetaEx(i), thetaEx(i)], [0, max(f)], 'r--')    
end