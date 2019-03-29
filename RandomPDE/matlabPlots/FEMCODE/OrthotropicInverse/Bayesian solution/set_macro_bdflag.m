function mesh = set_macro_bdflag(mesh)

neufun = @(x) (x(:,2) > mesh.box(4) - mesh.mesh_size*1e-3);
dirifun = @(x) (x(:,2) < mesh.box(1) + mesh.mesh_size*1e-3);
mesh.bdflag = cell(2,1);
mesh.bdflag{1,1} = get_bdflag(mesh, dirifun);
mesh.bdflag{2,1} = get_bdflag(mesh, neufun);

end