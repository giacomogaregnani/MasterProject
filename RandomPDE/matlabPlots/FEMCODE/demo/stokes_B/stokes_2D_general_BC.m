clear vp sol;
vp.a = 1;
vp.elemtype = 'p3';
vp.pelemtype = 'p2';

vp.f = @fstokesexB;
vp.bc = cell(2,2);

vp.bc{1,1} = @(x)(ex_stokesB(x)*[1;0]);
vp.bc{1,2} = @(x)(ex_stokesB(x)*[0;1]);
vp.bc{2,1} = @(x,normals)((ex_stokesB(x,[1,0])*[1;0]).*normals(:,1)+ ...
    (ex_stokesB(x,[0,1])*[1;0]).*normals(:,2)-ex_stokespB(x).*normals(:,1));
vp.bc{2,2} = @(x,normals)((ex_stokesB(x,[1,0])*[0;1]).*normals(:,1)+...
    (ex_stokesB(x,[0,1])*[0;1]).*normals(:,2)-ex_stokespB(x).*normals(:,2));

mesh = structured_mesh([0,1,0,1], 2, struct('centre',true));

neufun = @(x)(x(:,1) > 1-1e-5);
dirfun = @(x)(x(:,1) <= 1-1e-5);
mesh.bdflag{1,1} = get_bdflag(mesh, dirfun);
mesh.bdflag{1,2} = mesh.bdflag{1,1};
mesh.bdflag{2,1} = get_bdflag(mesh, neufun);
mesh.bdflag{2,2} = mesh.bdflag{2,1};

Nmax = 3;
[h1,l2,h1s,l2p,dof] = deal(zeros(Nmax,1));

for i=1:Nmax
	mesh = uniformrefine(mesh);
	[usol, ufemspace, psol, pfemspace, mesh, vp] = stokes(mesh, vp);	
	[h1(i), l2(i), h1s(i)] = get_H1error_exact(mesh, ufemspace, usol, @ex_stokesB);
	l2p(i) = get_L2error_exact( mesh, pfemspace, psol, @ex_stokespB);
	dof(i) = ufemspace.ndof * 2 + pfemspace.ndof;
end

figure;
loglog(dof,h1,dof,l2,dof,h1s,dof,l2p);


