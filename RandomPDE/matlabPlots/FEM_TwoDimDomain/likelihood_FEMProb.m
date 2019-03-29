function likelihood = likelihood_FEMProb(theta, V, D, xObs, uObs, obsNoise, MESHDet, data)

% Generate diffusion from theta
KL_Mean = zeros(size(V(:, 1)));
data.diffusion = exp(KL(KL_Mean, V, D, theta));

% Perturb vertices
MESH = MESHDet;
idxInternal = setdiff(1:size(MESH.vertices, 2), [MESH.boundaries(1, :), MESH.boundaries(2, :)]);
MESH.vertices(:, idxInternal) = MESH.vertices(:, idxInternal) + max(MESH.h)^2 * (-0.5 + rand(size(MESH.vertices(:, idxInternal))));

% Compute FEM solution
U = Elliptic_Solver(2, MESH.elements, MESH.vertices, MESH.boundaries, 'P1', data, [], [], false);

% Evaluate on observation points
uPoints = evalFEM(MESHDet, U, xObs);

% Compute log-likelihood
likelihood = -0.5 / (obsNoise^2) * (uPoints - uObs)' * (uPoints - uObs);

