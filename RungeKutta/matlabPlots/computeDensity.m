function [W, l, V] = computeDensity(fileName, tol)

V = dlmread([fileName, '.txt']);
dim = size(V, 2) / 2;
radius = V(1, dim + 1 : end);
volume = prod(radius);

x = dlmread([fileName, 'Matrix.txt'], ' ', 2, 0);
TransMatrix = sparse(x(:, 1), x(:, 2), x(:, 3));
TransMatrix = reshapeTransMatrix(TransMatrix);

% figure
% spy(TransMatrix)
% title('Transition Matrix sparsity pattern')

opts.tol = tol;
[W, l] = eigs(TransMatrix, 1, 'lm', opts);
W = abs(W);
densIntegral = sum(W * volume);
W = W / densIntegral;