% Test if the variance depends on h
clc; clear; close all;

getResults = false;

if getResults
    results = dlmread(['likelihood_EE_30_11_2016_08_43_.txt']);
    nExpEE = size(results, 1);
    hEE = zeros(nExpEE, 1);
    variances = hEE;
    biasSqd = hEE;
    
    for i = 1 : nExpEE
        hEE(i) = results(i, 1);
        variances(i) = results(i, 2);
        biasSqd(i) = results(i, 3);
    end
    
    results = dlmread(['likelihood_RK_30_11_2016_08_51_.txt']);
    
    nExpRK = size(results, 1);
    hRK = zeros(nExpRK, 1);
    variancesRK = hRK;
    biasSqdRK = hRK;
    
    for i = 1 : nExpRK
        hRK(i) = results(i, 1);
        variancesRK(i) = results(i, 2);
        biasSqdRK(i) = results(i, 3);
    end
else
    % 300 repetitions, final time T = 1, nData = 10
    load('likelihoodFITZNAG')    
end

markers = 'o*+<';
figure
loglog(hEE, variances, 'o-', 'marker', markers(1))
hold on
loglog(hRK, variancesRK, 'o-', 'marker', markers(2))
loglog(hEE, biasSqd, 'o-', 'marker', markers(3))
loglog(hRK, biasSqdRK, 'o-', 'marker', markers(4))
loglog(hEE, 100 * hEE.^2, 'k--')
loglog(hRK, 100 * hRK.^8, 'k')
xlim([1e-5 0.2])
legend('varEE', 'varRK', 'biasEE', 'biasRK', 'Location', 'SW')

orderEE = -log2(variances(2 : end) ./ variances(1 : end-1));
