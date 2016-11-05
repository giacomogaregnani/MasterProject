% Read results

clc; clear; close all

% Import data
results = dlmread('../data/orders.txt');
size = 2;
nExp = length(results) / size - 1; 

h = 0.25 ./ (2 .^ [0 : 3]);

refSol = results(1 : size);

res = zeros(size, nExp);
err = zeros(1, nExp);
for i = 1 : nExp
    res(:, i) = results(i * size + 1 : (i + 1) * size);
    err(i) = norm(res(:, i) - refSol);
end

loglog(h, err, 'o-');
hold on
loglog(h, 10 * h.^2, 'k');
loglog(h, 10 * h.^3, 'k--');
loglog(h, 10 * h.^4, 'k-.');
