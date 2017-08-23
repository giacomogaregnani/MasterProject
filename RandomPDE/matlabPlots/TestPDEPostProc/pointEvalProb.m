function uPoint = pointEvalProb(u, x, mesh, p)

uPoint = zeros(length(x), 1);

P = -1 / 2 + rand(1, mesh.N+1);
mesh.x = mesh.x + mesh.h^p * P;

uPoint = interp1(mesh.x, u, x)';