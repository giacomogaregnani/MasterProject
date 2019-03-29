function vp = get_variational_problem(pde,parameters)

vp.elemtype = 'p1';
vp.a =@(x,k,l) feval(pde.tensor,x,k,l,parameters);
vp.f = pde.f;
vp.bc{1} = @(p) pde.g_D(p);
vp.bc{2} = @(p, normals) pde.g_N(p,normals);

end