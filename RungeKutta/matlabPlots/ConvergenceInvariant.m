%% Parameters
clc; close all; clear;

nExp = 10;
h = 0.01 ./ (2.^[0 : nExp - 1]);

%% Read and compute error

E = zeros(nExp, 1);

figure
hold on
Inv = dlmread('IMConvergence.txt');
E = abs(Inv(:, 1) - Inv(:, 2));

figure
loglog(h, E);
