clc; clear; close all;

h = 0.1;
p = 1.5;

f1 = @(x) 1 ./ (x.^2 * h^p) .* (1 + 1 / (h^p) * ((1 ./ x) - h));
f2 = @(x) 1 ./ (x.^2 * h^p) .* (1 - 1 / (h^p) * ((1 ./ x) - h));

f = @(x) f1(x) .* (x > 1 / h) .* (x < 1 / (h - h^p)) + ...
         f2(x) .* (x < 1 / h) .* (x > 1 / (h + h^p));
     
X = linspace(1 / (h + h^p), 1 / (h - h^p), 10000);
plot(X, f(X))
hold on
plot([1/h, 1/h], [0, max(f(X))])