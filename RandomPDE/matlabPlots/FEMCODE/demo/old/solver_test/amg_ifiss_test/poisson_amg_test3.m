function resvec = poisson_amg_test3
%POISSON_AMG_TEST3 Summary of this function goes here
%   Detailed explanation goes here

mesh = structured_mesh([0,1,0,1,0,1], 16);

clear vp options;
vp.f = @f1;
vp.elemtype = 'p1';
vp.a = 1;
options.eqn = 'all';

%% Solve
[sol, femspace, ~, ~, ~, eqn] = poisson(mesh, vp, options);

%% Plotting

grid_data = amg_grids_setup(eqn.A);
smoother_params = amg_smoother_params(grid_data);
smoother_data = amg_smoother_setup(grid_data, smoother_params);

function pr = prec(x)
pr = amg_v_cycle(x, grid_data, smoother_data);
end


[sol2, ~, ~, ~, resvec] = pcg(eqn.A, eqn.f, 1e-8, 100, @prec);


NT = size(mesh.elem,1);
get_H1error(mesh, femspace, sol, mesh, femspace, sol2,1:NT)

end

