function femspace = p2_dof(mesh)
%% DOF_P2 dof structure for P2 element.

dim = size(mesh.elem,2) - 1;
N  = size(mesh.node,1);
NT = size(mesh.elem,1);

femspace.ldof = (dim+1)*(dim+2)/2;

if (dim == 1)
	bdnode = auxstructure(mesh,'bdnode');
	femspace.elem2dof = [];
	femspace.bd_dof = bdnode;
	femspace.ndof = N + NT;
    femspace.lastbdnode = N;
else
	[elem2edge, isbdedge, bdnode] = ...
		auxstructure(mesh, 'elem2edge','isbdedge','bdnode');
    NE = numel(isbdedge);
	femspace.elem2dof = N + double(elem2edge);
	femspace.bd_dof   = [bdnode; N + find(isbdedge)];
	femspace.ndof = N + NE;
end

femspace.deg  = 2;
femspace.elemtype = 'p2';
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