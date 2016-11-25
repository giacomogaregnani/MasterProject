% Test if the variance depends on h
clc; clear; close all;

results = dlmread('varMCResults.txt');
resultsRK = dlmread('varMCResults_RK.txt');

hEM = results(:, 1);
hRK = resultsRK(:, 1);
varEM = results(:, 2);
errEM = results(:, 3);
varRK = resultsRK(:, 2);
errRK = resultsRK(:, 3);

markers = 'o*+<';
loglog(hEM, varEM, 'marker', markers(1))
hold on
loglog(hRK, varRK, 'marker', markers(2))
loglog(hEM, errEM, 'marker', markers(3))
loglog(hRK, errRK, 'marker', markers(4))
loglog(hEM, hEM.^2, 'k')
loglog(hRK, hRK.^8, '--k')
legend('VarianceEM', 'VarianceRK', 'ErrorEM', 'ErrorRK', 'slope 2', 'slope 8', 'Location', 'SW')
xlabel('h')
ylabel('Var MC')

orderEM = -log2(varEM(2 : end) ./ varEM(1 : end-1));
orderRK = -log2(varRK(2 : end) ./ varRK(1 : end-1));