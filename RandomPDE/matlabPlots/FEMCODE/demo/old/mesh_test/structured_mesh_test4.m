%% generate mesh
mesh = structured_mesh([0,1],10);

%% display mesh
clear options;
options.nodenum = 1;
options.elemnum = 1;

figure;
simpplot(mesh, options);
