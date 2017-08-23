function mcL = likelihood_ProbFEM(theta, xObs, uObs, obsNoise, mE, mesh, data, nMC)
% Compute MC estimator of probabilistic log-likelihood function

L = zeros(1, nMC);

for i = 1 : nMC
    meshProb = mesh;
    meshProb.xInt = mesh.xInt + mesh.h^2 * (-0.5 + rand(size(mesh.xInt)));
    
    kappa = buildField_Andrea(theta);
    u = pointEval(solveFwdProblem_Cont(meshProb, kappa, data.f, data.rBC), xObs, mesh);
    
    L(i) = -obsNoise^(-2) * (u - uObs)' * (u - uObs);
end

% Avoid overflow in the sum of exponentials
[maxL, maxI] = max(L);
L(maxI) = [];

mcL = maxL + log(1 + sum(exp(L - maxL))) - log(nMC);