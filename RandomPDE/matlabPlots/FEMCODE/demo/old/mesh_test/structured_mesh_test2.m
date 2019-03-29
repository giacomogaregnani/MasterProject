%% generate mesh
clear options;
options.bc = 'periodic';
mesh = structured_mesh([0,1,0,1],[8,8], options);

%% display mesh
clear options;
options.nodenum = 1;
options.elemnum = 1;

figure;
simpplot(mesh,options);