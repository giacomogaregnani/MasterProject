function u = solveFwdProblemProb_Cont(mesh, kappa, f, rightBC)

tempVec = [mesh.xMin, mesh.xInt, mesh.xMax];
midPoints = 0.5 * (tempVec(1:end-1) + tempVec(2:end));
midPoints = [midPoints, 0];
gridVec = [tempVec; midPoints];
gridVec = gridVec(:)';
gridVec = gridVec(1:end-1);

kappaVec = kappa(gridVec);

A = assembleMatrixProb_Cont(kappaVec, mesh);
F = assembleRHSProb(f, mesh, rightBC);

% u = [0; A \ F; 0];
u = [0; A \ F];

end