function mesh = set_micro_bdflag(mesh)

mesh.bdflag = cell(2,1);
mesh.bdflag{1,1} = zeros(0,2);
mesh.bdflag{2,1} = zeros(0,2);

end