clc; clear; close all;

x = [0, sort(rand(1, 9)), 1];

A = sparse(zeros(9));

A(1, 1) = 1 / (x(2) - x(1)) + 1 / (x(3) - x(2));
A(1, 2) = -1 / (x(3) - x(2));
for i = 2 : 8
    A(i, i) = 1 / (x(i+1) - x(i)) + 1 / (x(i+2) - x(i+1));
    A(i, i-1) = - 1 / (x(i+1) - x(i));
    A(i, i+1) = - 1 / (x(i+2) - x(i+1));
end
A(9, 9) = 1 / (x(end-1) - x(end-2)) + 1 / (x(end) - x(end-1));
A(9, 8) = - 1 / (x(end-1) - x(end-2));

F = ones(9, 1);
for i = 2 : 10
    F(i-1) = F(i-1) * ((x(i+1) - x(i)) / 2 + (x(i) - x(i-1)) / 2);
end

xPlot = linspace(0, 1, 10000);
plot(x, [0; A \ F; 0])
hold on
plot(xPlot, -1 / 2 * (xPlot.^2 - xPlot))