function  femspace = p0d_dof( mesh )
%P0D_DOF Summary of this function goes here
%   Detailed explanation goes here
dim = size(mesh.elem, 2) -1;
NT = size(mesh.elem,1);
%femspace.bd_dof = auxstructure(mesh, 'bdnode'); NEEDS CORRECTION
femspace.ldof = 1;
femspace.ndof = femspace.ldof * NT; 
femspace.deg = 0;
femspace.elemtype = 'p0d';
femspace.lastbdnode = 0;
end

