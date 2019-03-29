clc; clear; close all
addpath('resultsErrorProb');
%%

y = dlmread('output.txt');
yDet = dlmread('outputdet.txt');

plot(y(:, 1), y(:, 2))
hold on
plot(yDet(:, 1), yDet(:, 2))