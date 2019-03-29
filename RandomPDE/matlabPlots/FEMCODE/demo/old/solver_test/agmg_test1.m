%% Problem definition
mesh = structured_mesh([0,1,0,1,0,1], 30);

clear vp;
vp.f = @f1;
vp.elemtype = 'p2';
vp.a = 1;
vp.solver='agmg';

%% Solve
[sol, femspace] = poisson(mesh, vp);

