clc; clear; close all
%%

T = 2000;
N = T * 10;
h = T / N;
tVec = linspace(0, T, N+1);

sigma = 2;
A = 3;
sqrt2Sigma = sqrt(2 * sigma);

nExp = 10000;

dW = sqrt(h) * randn(nExp, N);
X = zeros(nExp, N+1);

for i = 1 : N
    X(:, i+1) = X(:, i) - h * A * X(:, i) + sqrt2Sigma * dW(:, i);
end

% Compute the martingale
M = sum(X(:, 1:end-1) .* dW, 2);
% Compute the quadratic variation
qM = sum(X(:, 1:end-1).^2, 2) * h;

ratio = M ./ qM;

%%

results = sqrt(T) * sort(ratio);
histogram(results, 'normalization', 'pdf', 'numBins', 50);
hold on
plot(results, normpdf(results, 0, sqrt(A/sigma)));

figure
results = 1 / sqrt(T) * sort(M);
histogram(results, 'normalization', 'pdf', 'numBins', 50);
hold on
plot(results, normpdf(results, 0, sqrt(sigma/A)));