function [resvec] = poisson_amg_test2
%% Problem definition
fg.epsilon = 0.01;
fg.handle = @cg2_locper_inc_1;
sector.centre = [sqrt(2),sqrt(2)];
sector.delta = [fg.epsilon, fg.epsilon];

%% generate boundary
bdmesh = get_micro_geometry(fg, sector);

%% set approx geometry parameters
clear options;
options.startN = 8;
options.numiter = 15;
options.int_approx = false;
options.periodic = true;

%% generate approx. geometry
mesh = approx_mesh(bdmesh, options);
for i=1:1
	mesh = uniformrefine(mesh);
end

clear options;
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