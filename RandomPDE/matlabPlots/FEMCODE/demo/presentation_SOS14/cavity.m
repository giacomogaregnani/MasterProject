mesh = structured_mesh([-1,1,0,1], [200,100], struct('centre',true));
mesh.bdflag = 'dirichlet';
fh = @(x)([zeros(size(x,1),1), -ones(size(x,1),1)]);
bc1 = @(x)([double(x(:,2)==1), zeros(size(x,1),1)]);
bc2 = 0;
bc = {bc1, bc2};
vp = struct('a', 1, 'f', fh, ...
'elemtype', 'p2', 'pelemtype', 'p1');
vp.bc = bc;
[usol, ufemspace, psol, pfemspace, mesh, vp] = stokes(mesh, vp);
simpplot_sol(mesh, sqrt(sum(usol.^2,2)))
tricontour2(mesh, sqrt(sum(usol.^2,2)), 10)