function L = likelihood_FEM(theta, xObs, uObs, obsNoise, mesh, data)

kappa = buildField_Andrea(theta);
u = pointEval(solveFwdProblem_Cont(mesh, kappa, data.f, data.rBC), xObs, mesh);

L = -0.5 * obsNoise^(-2) * (u - uObs)' * (u - uObs);