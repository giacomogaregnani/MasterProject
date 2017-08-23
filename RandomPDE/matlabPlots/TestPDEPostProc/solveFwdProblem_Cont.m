function u = solveFwdProblem_Cont(mesh, kappa, f, rightBC)

% Interpolate on the grid
intGrid = linspace(mesh.xMin, mesh.xMax, 2 * mesh.N + 1);
kappaVec = kappa(intGrid);

A = assembleMatrix_Cont(kappaVec, mesh);
F = assembleRHS(f, mesh, rightBC);

u = [0; A \ F];

end