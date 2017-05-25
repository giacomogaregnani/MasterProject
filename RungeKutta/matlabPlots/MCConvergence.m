clc; clear; close all;

addpath('ResultsMSE')
results = dlmread('MSE_RK4_p5.txt');

h = results(:, 1);
MCErr = results(:, 2);

order = log2(MCErr(1:end-1)./MCErr(2:end));
meanOrder = mean(order);

% loglog(h, MCErr, 'o-', 'LineWidth', 2)
% hold on
% loglog(h, h.^0.5, 'r--')
% loglog(h, h, 'k--')
% loglog(h, h.^1.5, 'r')
% loglog(h, h.^2, 'k')
% loglog(h, h.^2.5, 'r.-')
% loglog(h, h.^3, 'k.-')
% loglog(h, h.^3.5, 'r.--')
% loglog(h, h.^5, 'k.--')
% legend('error', 'slope 0.5', 'slope 1', 'slope 1.5', 'slope 2', 'slope 2.5', 'slope 3', 'slope 3.5', 'slope 4', 'Location', 'SW')