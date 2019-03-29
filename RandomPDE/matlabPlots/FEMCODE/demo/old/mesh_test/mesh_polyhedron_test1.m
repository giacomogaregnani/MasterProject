%% define polyhedral geometry
clear bdmesh;
bdmesh.node = [0,0,0; 1,0,0; 0,1,0; -1,-1,1; -1,-1,-1];
bdmesh.elem = [1,4,2; 2,4,3; 3,4,1; 2,5,1; 3,5,2; 1,5,3];

%% define mesh parameters
clear options;
options.h = 0.1;
options.maxit = 50;
options.fix2 = true;

%% mesh geometry
mesh = mesh_polyhedron(bdmesh, options);

%% plot mesh
figure;
simpplot(mesh);

%% quality histogram
figure;
sq=simplex_quality(mesh);
hist(sq,0.01:0.01:0.99);

%% low quality 3D plot
figure;
mesh2=mesh;
mesh2.elem(sq>0.2,:) = [];
simpplot(mesh2);
simpplot(bdmesh);