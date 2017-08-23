function L = likelihood_CorrFEM(theta, xObs, uObs, obsNoise, mE, vE, mesh, data)

kappa = buildField_Andrea(theta);
u = pointEval(solveFwdProblem_Cont(mesh, kappa, data.f, data.rBC), xObs, mesh);

L = -0.5 * (u - uObs + mE)' * ((obsNoise^2 * eye(length(xObs)) + vE) \ (u - uObs + mE));