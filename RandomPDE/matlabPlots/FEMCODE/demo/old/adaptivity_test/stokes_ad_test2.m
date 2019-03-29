%% generate mesh
dim = 2;
if dim == 3
	mesh = structured_mesh([0,1,0,1,0,1],3,struct('bc','periodic'));
	ind = [5,11,13,14,15,17,23]-1;
	ind = [6*ind+1, 6*ind+2, 6*ind+3, 6*ind+4, 6*ind+5, 6*ind+6];
	mesh.elem(ind,:) = [];
elseif dim == 2
	mesh = structured_mesh([0,1,0,1],4,struct('bc','periodic'));
	mesh.elem(19,3) = 17;
	mesh.elem(21,2) = 17;
	mesh.elem(22,1) = 17;
	mesh.node = [mesh.node; mesh.node(11,:)];
end

%% variational problem
clear vp options;
vp.f = @fstokes;
vp.elemtype = 'p2';
vp.pelemtype = 'p1';
vp.a = 1;
vp.bc = 'zero_dirichlet';
%vp.solver = 'iuzawapcg';
%vp.iuzawapcgopt.verbose = true;

options.verbose = true;
%% adaptive process
vp.adapt.maxit = 50;
vp.adapt.maxdof = 10000;
vp.adapt.minabserr = 0.00001;
vp.adapt.minrelerr = 0.1;
vp.adapt.verbose = true;
%vp.adapt.marking.method = 'refine';

[usol, ufemspace, psol, pfemspace, mesh, vp] = ...
	stokes_ad(mesh, vp, options);

%% FINAL REFINEMENT
%[mesh,~,~,eqn{NI}.father] = uniformrefine(eqn{NI}.mesh);
%[usol, ufemspace, psol, pfemspace, mesh, vp] = ...
%	stokes(mesh, vp, options);


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
