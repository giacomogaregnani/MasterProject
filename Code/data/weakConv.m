clear; clc; close all;

resultsEE = dlmread('WeakConvergence.txt');
resultsMP = dlmread('WeakConvergence_MP.txt');
resultsRK = dlmread('WeakConvergence_RK.txt');
hEE = resultsEE(:, 1);
errEE = resultsEE(:, 2);
hRK = resultsRK(:, 1);
errRK = resultsRK(:, 2);
hMP = resultsMP(:, 1);
errMP = resultsMP(:, 2);

loglog(hEE, errEE, 'o-')
hold on
loglog(hRK, errRK, '*-')
loglog(hMP, errMP, '<-')
loglog(hEE, hEE, 'k--')
loglog(hMP, hMP.^2 / 3, 'k-.')
loglog(hRK, hRK.^4, 'k')
xlim([1e-4 0.2])
xlabel('h')

legend('errEE', 'errMP', 'errRK', 'slope 1', 'slope 2', 'slope 4', 'Location', 'SW');
