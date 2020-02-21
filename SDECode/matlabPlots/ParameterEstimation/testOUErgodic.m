clc; clear; close all;
%%

T = 5000;
N = 100 * T;
h = T / N;
timeVec = linspace(0, T, N+1);

V = @(x) x.^2/2;
F = @(x) x;
% V = @(x) x.^4/4 - x.^2/2;
% F = @(x) x^3 - x;
sigma = 0.5;
alpha = 1;
g = sqrt(2*sigma);

dW = sqrt(h) * randn(1, N);
W = [0, cumsum(dW, 2)];

% EM
X = zeros(1, N+1);
for i = 1 : N
    X(i+1) = X(i) + h * (-alpha * F(X(i))) + g * dW(i);
end

XX = linspace(min(X), max(X), 1000);
% theoryPdf = normpdf(XX, 0, sqrt(sigma));
theoryPdf = exp(-1/sigma * V(XX));
theoryPdf = theoryPdf / trapz(XX, theoryPdf);

figure
plot(timeVec, X);

figure
histogram(X, 'normalization', 'pdf')
hold on
plot(XX, theoryPdf, 'r')