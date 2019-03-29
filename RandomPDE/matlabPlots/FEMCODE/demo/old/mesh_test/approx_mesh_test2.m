%% set porous geometry parameters
clear fg sector options;
fg.epsilon = 0.1;
fg.handle = @cg3d_1;
sector.centre = [0.13,0.06,0.12];
sector.delta = [fg.epsilon, fg.epsilon, fg.epsilon];

%% generate boundary
bdmesh = get_micro_geometry(fg, sector);

%% set approx geometry parameters
clear options;
options.startN = 4;
options.numiter = 5;
options.intersects.approx = true;
options.periodic = true;

%% generate approx. mesh
mesh = approx_mesh(bdmesh, options);

%% plot
clear opt;
opt.wire = true;
opt.facealpha =0.1;
figure;
simpplot(bdmesh, opt);
simpplot(mesh);