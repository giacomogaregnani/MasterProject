clc; clear; close all
%% 

Y = dlmread('test.txt');

figure
plot(Y(:, 1), Y(:, 2))