%% generate mesh
mesh = structured_mesh([0,1,0,1,0,1],3,struct('bc','periodic'));
ind = [5,11,13,14,15,17,23]-1;
ind = [6*ind+1, 6*ind+2, 6*ind+3, 6*ind+4, 6*ind+5, 6*ind+6];
mesh.elem(ind,:) = [];

%options.lc = 0.1;
%[bdmesh, bdgmsh] = cg3d_2([1,1,1],0);
%mesh = gmsh(bdgmsh, options);

%options.lc = 0.2;
%[bdmesh, bdgmsh] = cg2d_1([1,1],0);
%mesh = gmsh(bdgmsh, options);
%mesh = label(mesh);

%% variational problem
clear vp options;
vp.f = @fstokes;
vp.elemtype = 'p1b';
vp.pelemtype = 'p1';
vp.a = 1;
vp.solver = 'iuzawapcg';
vp.iuzawapcgopt.verbose = true;

options.verbose = true;
%% adaptive process
vp.adapt.maxit = 300;
vp.adapt.maxdof = 2000000;
vp.adapt.minerr = 0.00001;
vp.adapt.verbose = true;
%vp.adapt.marking.method = 'refine';

[~, ~, ~, ~, ~, vp, stats, eqn] = stokes_ad(mesh, vp, options);
NI=stats.numiter;

%% FINAL REFINEMENT
[mesh,~,~,eqn{NI}.father] = uniformrefine(eqn{NI}.mesh);
[usol, ufemspace, psol, pfemspace, mesh, vp] = ...
	stokes(mesh, vp, options);


% %% Convergence rates
% offset = 0;
% figure;
% loglog(...
%     Ndof(1:end-offset),  uerrL2(1:end-offset,1,1),     '-x', ...
%     Ndof(1:end-offset),  uerrsemiH1(1:end-offset,1,1), '-o', ...
%     Ndof(1:end-offset),  errest(1:end-offset,1,1), '-o', ...
%     Ndof(1:end-offset),  perrL2(1:end-offset,1,1),     '-o');
% hold on;
% params={'interpreter','latex','FontSize',14};
% legend({'$\| u-u^h\|_{L^2}$', ...
%     '$| u-u^h|_{H^1}$', ...
% 	'error est', ...
%     '$\| p-p^h\|_{L^2}$'},params{:});
% xlabel('$N_{dof}$',params{:})
