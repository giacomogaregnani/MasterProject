function likelihood = likelihood_FEM(theta, uMean, V, D, xObs, uObs, obsNoise, MESH, data)

% Generate diffusion from theta
data.diffusion = exp(KL(uMean, V, D, theta));

% Compute FEM solution
U = Elliptic_Solver(2, MESH.elements, MESH.vertices, MESH.boundaries, 'P1', data, [], [], false);

% Evaluate on observation points
uPoints = evalFEM(MESH, U, xObs);

% Compute log-likelihood
likelihood = -0.5 / (obsNoise^2) * (uPoints - uObs)' * (uPoints - uObs);
