clc; clear; close all
%% 

trueVal = 1 - normcdf(4, 0, 1);

phi = @(x) x > 4;

N = 1000;

ZMC = 0;
for i = 1 : N
    X = randn;
    ZMC = ZMC + phi(X) / N;
end

ZIS = 0;
for i = 1 : N
    X = 2 + randn;
    ZIS = ZIS + 1 / N * phi(X) * normpdf(X) / normpdf(X, 2, 1);
end