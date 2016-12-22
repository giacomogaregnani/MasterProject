
clc; close all; clear

h = 0.1 ./ 2.^[0 : 5];
nExp = length(h);

errP = zeros(size(h));
errD = zeros(size(h));

trueVal = [0.2, 0.2, 3.0]';
phiTrueVal = trueVal' * trueVal;

for i = 1 : nExp
    str = ['allThetas_EE_30_11_2016_11_02_', num2str(h(i)*1e6)];
    resultsP = dlmread([str, '.txt']);
    meanP = mean(resultsP)';
    errP(i) = abs(meanP'*meanP - phiTrueVal);
    resultsD = dlmread([str, '_det.txt']);
    meanD = mean(resultsD)';
    errD(i) = abs(meanD'*meanD - phiTrueVal);
%     plotMultiRV(resultsP, trueVal, resultsD);
end

figure
loglog(h, errP, 'o-')
hold on
loglog(h, errD, 'o-')
loglog(h, 100 * h, 'k-')
