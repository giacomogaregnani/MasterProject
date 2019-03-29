%% define polyhedral geometry
clear bdmesh;
bdmesh.node = [0,0; 0,2; 1,2; 0.7,0.7; 2,1; 2,0];
bdmesh.elem = [1,6; 6,5; 5,4; 4,3; 3,2; 2,1];

%% define mesh parameters
clear options;
options.h = 0.2;
options.maxit = 55;
options.fix1 = true;

%% mesh geometry
mesh = mesh_polyhedron(bdmesh, options);

%% plot mesh
figure;
simpplot(mesh);

%% quality histogram
figure;
hist(simplex_quality(mesh),0.01:0.01:0.99);