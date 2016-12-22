
clear; close all; clc

nMCMC = [1000, 2000, 4000, 8000, 16000, 32000, 64000];
nExp = length(nMCMC);

for i = 1 : nExp
    result = dlmread(['SigmaExp/test_22_12_2016_09_58_', num2str(nMCMC(i)), '.txt']);
    means = result(:, 2);
    directVar(i) = var(means);
    sigmasBM = result(:, 1);
    meanSigmaBM(i) = mean(sigmasBM);
end
meanSigmaBM = meanSigmaBM ./ nMCMC;


figure
loglog(nMCMC, directVar, 'o-')
hold on
loglog(nMCMC, meanSigmaBM, '*-')
loglog(nMCMC, nMCMC.^(-1), 'k--')
legend('var', 'BM', 'slope -1')
xlabel('N')

figure
loglog(nMCMC, abs(directVar - meanSigmaBM), 'o-')
hold on
loglog(nMCMC, 100 * nMCMC.^-2, 'k--')
xlabel('N')
