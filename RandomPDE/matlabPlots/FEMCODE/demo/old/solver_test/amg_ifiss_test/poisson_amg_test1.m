function [eqn] = poisson_amg_test1

%% Problem definition
mesh = structured_mesh([0,1,0,1],128);

clear vp options;
vp.f = @f1;
vp.elemtype = 'p1';
vp.a = [2, 1; 0.2, 1];
options.eqn = 'all';

%% Solution
[sol, femspace, ~, ~, ~, eqn] = poisson(mesh, vp, options);

%% Plotting
simpplot_sol(mesh,sol);

grid_data = amg_grids_setup(eqn.A);
smoother_params = amg_smoother_params(grid_data);
smoother_data = amg_smoother_setup(grid_data, smoother_params);

function pr = prec(x)
pr = amg_v_cycle(x, grid_data, smoother_data);
end


[sol2, ~, ~, ~, resvec] = pcg(eqn.A, eqn.f, 1e-8, 100, @prec);


NT = size(mesh.elem,1);
get_H1error(mesh, femspace, sol, mesh, femspace, sol2,1:NT)

% for i=1:10
% 	get_H1error(mesh, femspace, sol, mesh, femspace, sol_amg,1:NT)
% 	sol_amg = amg_v_cycle(eqn.rhs, grid_data, smoother_data, i, sol_amg);
% end
end