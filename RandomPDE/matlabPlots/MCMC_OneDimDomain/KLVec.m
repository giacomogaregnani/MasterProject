clc; clear; close all

V = dlmread('TestKL.txt');

X = linspace(0, 1, 100);

for i = 1 : 5
    figure
    plot(X, V(i, :))
end