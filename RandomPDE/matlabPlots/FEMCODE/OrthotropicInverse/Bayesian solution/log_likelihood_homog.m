function l = log_likelihood_homog(param, mesh, locations, observations, std_noise)

pde = define_pde();
vp = get_variational_problem(pde, [exp(param(1)), param(2)]);
[sol, femspace, mesh] = my_elasticity_2d(mesh, vp);

results = get_observations(mesh, sol, locations, femspace);

d = results(:) - observations;

if size(std_noise, 1) == 1
    l = -1 / (2 * std_noise^2) * d'*d;
else
    l = -1 / 2 * d' * (std_noise \ d);
end

return