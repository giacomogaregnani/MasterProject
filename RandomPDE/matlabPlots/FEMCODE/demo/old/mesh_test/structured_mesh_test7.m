%% generate mesh
mesh = structured_mesh([0,1,0,1,0,1],1);

%% display mesh
clear opt;
opt.wire = true;
opt.facealpha = 0;
opt.nodenum = true;
opt.elemnum = true;
opt.fontsize = 25;
opt.facealpha = 0.1;
opt.view = 3;
figure;
simpplot(mesh, opt);


