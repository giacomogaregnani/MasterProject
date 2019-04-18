clc; clear; close all
%%

% obs = dlmread('testSol1.txt');
% PF = dlmread('test1.txt');
% PF2 = dlmread('test1.txt');
% 
% t = linspace(0, 10, length(obs));
% 
% figure
% plot(t, PF(1:end-1,:), 'color', 0.6 * ones(1, 3))
% hold on
% plot(t, PF(end,:), 'r')
% plot(t, obs, 'k')
% 
% figure
% plot(t, PF2(1:end-1,:), 'color', 0.6 * ones(1, 3))
% hold on
% plot(t, PF2(end,:), 'r')
% plot(t, obs, 'k')

l = dlmread('testLik.txt');
l2 = dlmread('testLik2.txt');
figure
hold on
H = histogram(l, 'normalization', 'pdf');
H2 = histogram(l2, 'normalization', 'pdf');
hold on
meanl = mean(l);
plot([meanl, meanl], [0, max(H.Values)], 'k--', 'linewidth', 2)
meanl2 = mean(l2);
plot([meanl2, meanl2], [0, max(H2.Values)], 'k', 'linewidth', 2)