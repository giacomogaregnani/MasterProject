% Test if the variance depends on h
clc; clear; close all;

resultsLoad = true;
MorH = 'h';

if resultsLoad
    results = dlmread('varMCResults.txt');
    resultsRK = dlmread('varMCResults_RK.txt');
    resultsBE = dlmread('varMCResults_BE_Lorenz.txt');
    
    hEM = results(:, 1);
    hRK = resultsRK(:, 1);
    hBE = resultsBE(:, 1);
    varEM = results(:, 2);
    errEM = results(:, 3);
    varRK = resultsRK(:, 2);
    errRK = resultsRK(:, 3);
    varBE = resultsBE(:, 2) / 800;
    errBE = resultsBE(:, 3) / 800;
else
    load('MCVarianceFITZNAG');
end

markers = 'o*+<>.';
if MorH == 'h'
    loglog(hEM, varEM, 'marker', markers(1))
    hold on
    loglog(hRK, varRK, 'marker', markers(2))
    loglog(hBE, varBE, 'marker', markers(6))
    loglog(hEM, errEM, 'marker', markers(3))
    loglog(hRK, errRK, 'marker', markers(4))
    loglog(hBE, errBE, 'marker', markers(5))
    loglog(hBE, 100000000000000 * hBE.^8, 'k')
    loglog(hRK, hRK.^8, '--k')
    legend('VarianceEM', 'VarianceRK', 'VarianceBE', 'ErrorEM', 'ErrorRK', 'ErrorBE', 'slope 2', 'slope 8', 'Location', 'SW')
    xlabel('h')
    ylabel('Var MC')
    
    orderEM_V = computeOrder(hEM, varEM);
    orderRK_V = computeOrder(hRK, varRK);
    orderBE_V = computeOrder(hBE, varBE);
    orderEM_B = computeOrder(hEM, errEM);
    orderRK_B = computeOrder(hRK, errRK);
    orderBE_B = computeOrder(hBE, errBE);
else
    markers = 'o*+<';
    loglog(1./hEM, varEM, 'marker', markers(1))
    hold on
    loglog(1./hRK, varRK, 'marker', markers(2))
    loglog(1./hEM, 1./hEM, 'k')
    legend('VarianceEM', 'VarianceRK', 'VarianceBE', 'slope 1', 'Location', 'SE')
    xlabel('M^{-1}')
    
    orderEM = -log2(varEM(2 : end) ./ varEM(1 : end-1));
    orderRK = -log2(varRK(2 : end) ./ varRK(1 : end-1));
end
