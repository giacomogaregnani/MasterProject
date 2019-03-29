function [V, D] = laplEig(x, nEig)

N = length(x);
L = x(end) - x(1);

V = zeros(N, nEig);
D = zeros(nEig, 1);

for i = 1 : nEig
   V(:, i) = sin(i * pi * x / L)';
   D(i) = -i^2 * pi^2;
end
V = V * sqrt(2 / L);
D = D / L^2;