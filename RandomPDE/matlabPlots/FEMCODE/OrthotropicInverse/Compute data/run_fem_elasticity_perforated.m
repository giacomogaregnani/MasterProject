%% Compute  full fine scale solution
clear
clc
tic
%% Define pde

%% Geometry and mesh
box = [0, 1, 0, 1]; theta = 0.8; epsilon = 50; radius = sqrt((1-theta)/pi)/epsilon;
mesh_size = 2*pi*radius/20;
pde = define_pde();
numObservations_perEdge = 10;
nr_edges = 3;

Observations = zeros(1, nr_edges, numObservations_perEdge, 2);
locations = linspace(box(1), box(2), numObservations_perEdge + 2);
locations = locations(2:end-1);
E = 7; nu = 0.3;
options = struct('dim', 2, 'verbose', 0);

figure
hold on
N = 500;

mesh = structured_mesh([0, 1, 0, 1], N, []);
mesh.box = box; mesh.mesh_size = mesh_size;
mesh = set_macro_bdflag(mesh);

%% Define variational problem
vp = get_variational_problem(pde,[E, nu]);
[sol, femspace, mesh, vp] = my_elasticity_2d(mesh, vp);
Observations(1,:,:,:) = get_observations(mesh, sol, locations, femspace);
plot(locations, reshape(Observations(1,2,:,2), numObservations_perEdge, 1), '-o')

filename = ['observations_numObs', num2str(numObservations_perEdge),'.mat'];
save(filename);

toc