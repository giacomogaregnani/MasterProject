clear; clc; close all;

results = dlmread('WeakConvergence.txt');
h = results(:, 1);
err = results(:, 2);

loglog(h, err, 'o-')
hold on
loglog(h, h, 'k--')