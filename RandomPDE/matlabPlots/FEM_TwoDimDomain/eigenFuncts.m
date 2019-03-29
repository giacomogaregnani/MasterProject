function [V, D] = eigenFuncts(MESH, NKL)

gamma = 0.05;
lambda = 0.5;

[a1, b1] = meshgrid(MESH.nodes(1, :)); [a2, b2] = meshgrid(MESH.nodes(2, :));
L = gamma * exp(-sqrt((a1 - b1).^2 + (a2 - b2).^2)/lambda);

[V, D] = eigs(L, NKL, 'lm');

end