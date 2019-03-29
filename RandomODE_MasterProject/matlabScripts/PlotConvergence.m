clc; clear; close all
%% 

data = dlmread('outputWeak.txt');

h = data(:, 1); 
err = data(:, 2);

loglog(h, err, 'ko-', 'linewidth', 1.5)
hold on
loglog(h, h.^(0.5), 'k--.')
loglog(h, h, 'k')
loglog(h, h.^(1.5), 'k.-')
loglog(h, h.^2, 'k--')
legend('err', 'h^{0.5}', 'h', 'h^{1.5}', 'h^{2}', 'location', 'best')