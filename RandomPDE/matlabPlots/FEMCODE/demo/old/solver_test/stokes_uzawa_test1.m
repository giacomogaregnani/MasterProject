%% set porous geometry parameters
load('mesh2d_stokes_per_small.mat')

%% set variational problem
clear vp;
vp.f = @fstokesmicro;
vp.elemtype = 'p2';
vp.pelemtype = 'p1';
vp.a = 1;
vp.solver = 'uzawa';

%% Solve
[usol, ufemspace, psol, pfemspace] = stokes(mesh, vp);

%% Plot pressure
report_stokes(mesh, vp, ufemspace, usol, pfemspace, psol);