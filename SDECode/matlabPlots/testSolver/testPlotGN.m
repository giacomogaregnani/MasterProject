clc; clear; close all
%% 

Y = dlmread('testGN.txt');

Y1 = [Y(1:101, 1), Y(1:101, 2)];
Y2 = [Y(102:1102, 1), Y(102:1102, 2)];
Y3 = [Y(1103:end, 1), Y(1103:end, 2)];

figure
plot(Y1(:, 1), Y1(:, 2))
hold on
plot(Y2(:, 1), Y2(:, 2))
plot(Y3(:, 1), Y3(:, 2))

Y = dlmread('test.txt');

plot(Y(:, 1), Y(:, 2), 'k')