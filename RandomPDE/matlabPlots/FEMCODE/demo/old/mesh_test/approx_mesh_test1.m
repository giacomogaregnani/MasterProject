%% set porous geometry parameters
clear fg sector options;
fg.epsilon = 0.1;
fg.handle = @cg2d_1;
sector.centre = [0.47,0.22];
sector.delta = fg.epsilon;

%% generate boundary
bdmesh = get_micro_geometry(fg, sector);

%% set approx geometry parameters
clear options;
options.startN = 8;
options.numiter = 1;
options.int_approx = false;
options.periodic = true;

%% generate approx. geometry
mesh = approx_mesh(bdmesh, options);

%% plot
figure;
simpplot(mesh);
simpplot(bdmesh);
