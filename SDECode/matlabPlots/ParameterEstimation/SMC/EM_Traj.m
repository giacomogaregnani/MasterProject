function [t, Y] = EM_Traj(y0, f, g, T, N, M)
% [t, y] = EM_UO(y0, a, S, T, N, M) 
% Approximates the solution of dY = f(Y) dt + g(Y) dWt
% with Euler-Maruyama

h = T / N;
Y = zeros(M, N+1);
t = linspace(0, T, N+1);
Y(:, 1) = y0;

for j = 1 : N
    Y(:, j+1) = EM_1S(Y(:, j), f, g, h);
%     xi = randn(M, 1); 
%     Y(:, j+1) = Y(:, j) + h * f(Y(:, j)) + sqrt(h) * g(Y(:, j)) .* xi;
end
