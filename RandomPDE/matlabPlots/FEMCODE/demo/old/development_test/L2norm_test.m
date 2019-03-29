%% Problem definition
mesh = structured_mesh([0,1,0,1],32);

clear vp;
vp.f = @f1;
vp.a = 1;

%% Solution
vp.elemtype = 'p1b';
[sol, femspace] = poisson(mesh, vp);
get_L2norm(mesh, femspace, sol)

vp.elemtype = 'p1';
[sol, femspace] = poisson(mesh, vp);
get_L2norm(mesh, femspace, sol)

vp.elemtype = 'p2';
[sol, femspace] = poisson(mesh, vp);
get_L2norm(mesh, femspace, sol)