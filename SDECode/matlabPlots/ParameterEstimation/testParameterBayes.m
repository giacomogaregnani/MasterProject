clc; clear; close all
%% 

results = dlmread('driftBayesian.txt');

%%
truth = results(1);
MLE = results(2);
post = results(3:end);

[f, xi] = ksdensity(post);

plot(xi, f, 'k-')
hold on
plot([truth, truth], [0, max(f)], 'k--')
plot([MLE, MLE], [0, max(f)], 'k-.')
legend('post', 'truth', 'MLE')