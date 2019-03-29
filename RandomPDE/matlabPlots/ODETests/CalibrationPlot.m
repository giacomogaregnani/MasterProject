clc; clear; close all
addpath('resultsCalibration');
%%

filename = 'calibPendulumHam';

figure 
hold on

for i = 1 : 1
    x = dlmread([filename '.txt']);
    
    [xSort, I] = sort(x(:, 1));
    densities = exp(x(I, 2)) / trapz(xSort, exp(x(I, 2)));
    plot(xSort, densities);
    hold on
    [f, xi] = ksdensity(x(:, 1));
    plot(xi, f, 'linewidth', 1)
end

figure
y = dlmread([filename 'Growth.txt']);
N = size(y, 1);
loglog(y(1:N/2, 1), y(1:N/2, 2))
hold on
loglog(y(N/2+1:N, 1), y(N/2+1:N, 2))

y = dlmread([filename 'Growth2.txt']);
loglog(y(:, 1), y(:, 2))
loglog(y(2:end, 1), abs(y(2:end, 3) - y(1, 3)))

loglog(y(:, 1), 0.5e-4*y(:, 1).^2, 'k-.')
loglog(y(:, 1), 1e-4*y(:, 1), 'k--')
loglog(y(:, 1), 2e-4*y(:, 1).^(1/2), 'k')
loglog(y(:, 1), 0.5e-4*y(:, 1).^(3/2), 'k-.')

legend('error', 'stddev1', 'stddev2', 'location', 'SE')