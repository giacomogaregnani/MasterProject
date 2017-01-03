
clear; close all; clc

hText = [20000, 10000, 5000, 2500, 1250, 625];
h = hText * 1e-5;
nExp = length(hText);

Z = norminv(0.975);
nBatchTextes = floor(sqrt(50000));
trueVal = norm([0.2, 0.2, 3.0])^2;

for i = 1 : nExp
    result = dlmread(['ProbVsDet/testOne_02_01_2017_11_13_', num2str(hText(i)), 'prob.txt']);
    means(i) = result(2);
    sigmaBMProb(i) = result(1);
end
quantileProb = Z * sigmaBMProb / sqrt(nBatchTextes);
intervalProb = [means' - quantileProb', means' + quantileProb'];

errorbar(h, means, quantileProb)

for i = 1 : nExp
    result = dlmread(['ProbVsDet/testOne_02_01_2017_11_13_', num2str(hText(i)), 'det.txt']);
    means(i) = result(2);
    sigmaBM(i) = result(1);
end
quantileDet = Z * sigmaBM / sqrt(nBatchTextes);
intervalDet = [means' - quantileDet', means' + quantileDet'];

hold on
errorbar(h, means, quantileDet)

% figure
% loglog(h, intervalDet)
% hold on
% loglog(h, intervalProb)
% loglog(h, repmat(trueVal, nExp, 1))
