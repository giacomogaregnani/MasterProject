function l = log_likelihood_homog_prob_2(param, mesh, locations, observations, std_noise, nMC)

pde = define_pde();
vp = get_variational_problem(pde, [exp(param(1)), param(2)]);

N_DOF = size(mesh.node, 1);
N = sqrt(N_DOF);
[sol, femspace, ~] = my_elasticity_2d(mesh, vp);

l = zeros(nMC, 1);
for i = 1 : nMC
    solProb = sol + 1 / N * randn(size(sol));
    results = get_observations(mesh, sol, locations, femspace);
    d = results(:) - observations(:);
    l(i) = -1 / (2 * std_noise^2) * d'*d;
end

maxL = max(l);
l = l(setdiff(1:nMC, maxL));
sumL = sum(exp(l - maxL));

l = maxL + log(1 + sumL) - log(nMC);


return