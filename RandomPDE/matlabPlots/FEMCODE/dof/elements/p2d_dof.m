function  femspace = p2d_dof(mesh)
%P1D_DOF Summary of this function goes here
%   Detailed explanation goes here

dim = size(mesh.elem, 2) -1;
NT = size(mesh.elem,1);

femspace.ldof = (dim+1) + (dim + 1)* dim/2;
femspace.ndof = femspace.ldof * NT;
femspace.deg = 2;
femspace.elemtype = 'p2d';

if (dim == 1)
	femspace.nodelambda = [eye(dim+1); 1/2*[1,1]];
elseif (dim == 2)
	femspace.nodelambda = [eye(dim+1); 1/2 * [1,1,0; 1,0,1; 0,1,1]];
elseif (dim == 3)
	femspace.nodelambda = [eye(dim+1); 1/2 * ...
		[1,1,0,0; 1,0,1,0; 1,0,0,1; 0,1,1,0; 0,1,0,1; 0,0,1,1]];
else
	error('nodelambda not defined')
end

end