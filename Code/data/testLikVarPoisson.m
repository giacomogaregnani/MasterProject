% Test if the variance depends on h
clc; clear; close all;

hEE = 0.1 ./ 2.^[0:12];
hTextEE = floor(hEE * 1e9);

nExpEE = length(hEE);
variances = zeros(size(hEE));
bias = zeros(size(hEE));

exact =  dlmread(['likelihoodVar_poisson_24_11_2016_01_49_exact.txt']);

for i = 1 : nExpEE
    results = dlmread(['likelihoodVar_poisson_24_11_2016_01_49__h_', num2str(hTextEE(i)), '.txt']);
    variances(i) = var(results);
    bias(i) = abs(mean(results) - exact)^2;
end

figure
loglog(hEE, variances, 'o-')
hold on
loglog(hEE, bias, 'o-')
loglog(hEE, 1000 * hEE.^2, 'k--')
loglog(hEE, 1000 * hEE.^4, 'k')

orderEE = -log2(variances(2 : end) ./ variances(1 : end-1));
orderEEbias =-log2(bias(2 : end) ./ bias(1 : end-1));