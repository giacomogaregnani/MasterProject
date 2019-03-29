function  femspace = p1d_dof( mesh )
%P1D_DOF Summary of this function goes here
%   Detailed explanation goes here

dim = size(mesh.elem, 2) -1;
NT = size(mesh.elem,1);

femspace.ldof = dim + 1;
femspace.ndof = femspace.ldof * NT; 
femspace.deg = 1;
femspace.elemtype = 'p1d';

femspace.nodelambda = [eye(dim+1); repmat(1/(dim+1),[1,dim+1])];
