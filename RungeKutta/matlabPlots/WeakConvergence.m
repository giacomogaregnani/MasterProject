clc; clear; close all;

SolutionRandH = dlmread('weakConv.txt');
h = 0.1 ./ 2 .^ [0 : length(SolutionRandH) - 1];

% nMC = 10000;
% C = 1 / 3;

% p = 1.5;

% for i = 1 : length(SolutionRandH)
% 
%     Brownian = randn(1, 10000);
%     Y = exp(1 - 0.5 * C * h(i)^(2 * p - 1) + sqrt(C * h(i)^(2 * p - 1)) * Brownian);
%     MeanY(i) = mean(Y);
%     
%     err(i) = abs(SolutionRandH(i, 1) - MeanY(i));
% end

weakErr = abs(SolutionRandH(:, 1) - SolutionRandH(:, 2));

order = log2(weakErr(1:end-1)./weakErr(2:end));
meanOrder = mean(order);

loglog(h, weakErr, 'o-')
hold on
loglog(h, h, 'k--')
loglog(h, h.^2, 'k')
loglog(h, h.^3, 'k.-')
loglog(h, h.^4, 'k.--')
legend('error', 'slope 1', 'slope 2', 'slope 3', 'slope 4', 'Location', 'SE')
