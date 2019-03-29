function l = log_likelihood_homog_prob(param, mesh, locations, obsVec, std_noise, nMC)

pde = define_pde();
vp = get_variational_problem(pde, [exp(param(1)), param(2)]);

N_DOF = size(mesh.node, 1);
N = sqrt(N_DOF);

femspace = get_femspace(mesh, 'p1');
idxInt = setdiff(1:N_DOF, femspace.bd_dof);
intNodes = mesh.node(idxInt, :);
nInt = size(intNodes, 1);

l = zeros(nMC, 1);

for i = 1 : nMC
    MESH = mesh;
    pertNodes = intNodes + 0.5 / N * (-0.5 + rand(nInt, 2));
    MESH.node(idxInt, :) = pertNodes;
    
    [sol, femspace, ~] = my_elasticity_2d(MESH, vp);
    results = get_observations(mesh, sol, locations, femspace);
    
    d = results(:) - obsVec;
    l(i) = -1 / (2 * std_noise^2) * (d'* d);
end

maxL = max(l);
l = l(setdiff(1:nMC, maxL));
sumL = sum(exp(l - maxL));

l = maxL + log(1 + sumL) - log(nMC);


return