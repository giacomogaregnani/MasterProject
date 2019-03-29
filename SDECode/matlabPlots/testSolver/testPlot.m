clc; clear; close all;
%%

X = dlmread('test.txt');

x = linspace(-3, 3, 1000);
V = x.^4 / 4 - x.^2 / 2;
rho = exp(-V / 0.5); 
rho = rho / trapz(x, rho);

histogram(X, 'normalization', 'pdf')
hold on
plot(x, rho, 'k', 'linewidth', 2)