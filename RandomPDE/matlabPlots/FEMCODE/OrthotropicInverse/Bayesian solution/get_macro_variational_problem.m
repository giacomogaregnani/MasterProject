function vp = get_macro_variational_problem(pde)

vp.elemtype = 'p1';
vp.f = pde.f;
vp.bc{1} = @(p) pde.g_D(p);
vp.bc{2} = @(p, normals) pde.g_N(p,normals);

end