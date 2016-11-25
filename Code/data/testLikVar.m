% Test if the variance depends on h
clc; clear; close all;

hRK = 0.1 ./ 2.^[0:9];
hEE = 0.1 ./ 2.^[0:12];
hTextRK = floor(hRK * 1e9);
hTextEE = floor(hEE * 1e9);

nExpRK = length(hRK);
nExpEE = length(hEE);
variances = zeros(size(hEE));
variancesRK = zeros(size(hRK));

for i = 1 : nExpRK
    resultsRK = dlmread(['likelihood_var_RK_24_11_2016_09_58__h_', num2str(hTextRK(i)), '.txt']);
    variancesRK(i) = var(resultsRK);
end

for i = 1 : nExpEE
    results = dlmread(['likelihood_var23_11_2016_05_46__h_', num2str(hTextEE(i)), '.txt']);
    variances(i) = var(results(:, 4));
end

figure
loglog(hEE, variances, 'o-')
hold on
loglog(hRK, variancesRK, 'o-')
loglog(hEE, 1000 * hEE.^2, 'k--')
loglog(hRK, 1000 * hRK.^8, 'k')

orderEE = -log2(variances(2 : end) ./ variances(1 : end-1));
orderRK = -log2(variancesRK(2 : end) ./ variancesRK(1 : end-1));