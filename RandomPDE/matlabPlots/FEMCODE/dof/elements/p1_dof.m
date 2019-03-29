function femspace = p1_dof(mesh)
%P1_DOF Summary of this function goes here
%   Detailed explanation goes here

dim = size(mesh.elem, 2) -1;

femspace.bd_dof = auxstructure(mesh, 'bdnode');
femspace.ndof = size(mesh.node,1);
femspace.ldof = dim + 1;
femspace.deg = 1;
femspace.elemtype = 'p1';
femspace.nodelambda = eye(dim+1);
end

