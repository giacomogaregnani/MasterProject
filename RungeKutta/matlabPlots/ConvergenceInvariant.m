clc; close all; clear;

%% Read and compute error

results = dlmread('invariantError.txt');

loglog(results(:, 1), results(:, 2), 'o-')
hold on
loglog(results(:, 1), results(:, 1).^2, 'k--')