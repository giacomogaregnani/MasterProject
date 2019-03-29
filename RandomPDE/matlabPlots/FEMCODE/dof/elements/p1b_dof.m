function femspace = p1b_dof(mesh)
%P1B_DOF2 Summary of this function goes here
%% DOF_P2 dof structure for P2 element.

N = size(mesh.node,1);
NT = size(mesh.elem,1);
dim = size(mesh.elem, 2) -1;

femspace.elem2dof = (N+1 : N+NT)';

femspace.bd_dof = auxstructure(mesh, 'bdnode');
femspace.ndof = N+NT;
femspace.ldof = dim + 2;
femspace.deg = dim+1;
femspace.elemtype = 'p1b';
femspace.nodelambda = [eye(dim+1); repmat(1/(dim+1),[1,dim+1])];