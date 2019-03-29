clc; clear; close all
%%

nKL = 10;

N = 20;
x = linspace(0, 1, N);
y = linspace(0, 1, N);
[XX, YY] = meshgrid(x, y);

gamma = 0.05;
lambda = 0.01;

[a1, b1] = meshgrid(repmat(x, 1, length(y))); 
[a2, b2] = meshgrid(reshape(repmat(y, length(x), 1), 1, N^2));
L = gamma * exp(-sqrt((a1 - b1).^2 + (a2 - b2).^2)/lambda);

[V, D] = eigs(L, nKL, 'lm');

for i = 1 : nKL
    figure
    eFunc = reshape(V(:, i), N, N);  
    surf(XX, YY, eFunc, 'edgecolor', 'none')
    colormap('jet')
    view(2)
end
