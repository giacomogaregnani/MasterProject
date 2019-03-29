clc; clear; close all;
%%

x = dlmread('SMC1.txt');
x2 = dlmread('SMC21.txt');
y = dlmread('MCMC_ODE21.txt');
thetaEx = [0.2, 0.2, 3.0];

for i = 1 : length(thetaEx)
    figure
    [f, xi] = ksdensity(exp(x(:, i+1)), 'weights', x(:, 1));
    plot(xi, f)  
    hold on
    [f, xi] = ksdensity(exp(x2(:, i+1)), 'weights', x2(:, 1));
    plot(xi, f)  
    [f, xi] = ksdensity(exp(y(:, i)));
    plot(xi, f)    
    plot(([thetaEx(i), thetaEx(i)]), [0, max(f)], 'k--')    
    legend('prob', 'det', 'MCMC', 'truth')
end