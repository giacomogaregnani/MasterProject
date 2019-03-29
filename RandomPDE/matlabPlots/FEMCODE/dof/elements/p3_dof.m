function femspace = p3_dof(mesh)
%% DOF_P2 dof structure for P2 element.

dim = size(mesh.elem,2) - 1;
N  = size(mesh.node,1);
NT = size(mesh.elem,1);

femspace.ldof = ...
    (dim+1) + ...
    2 * (dim+1)*dim/2 + ...
    1 * (dim+1)*dim*(dim-1)/6;

if (dim == 1)
	bdnode = auxstructure(mesh,'bdnode');
	femspace.elem2dof = [];
	femspace.bd_dof = bdnode;
	femspace.ndof = N + 2*NT;
    femspace.lastbdnode = N;
elseif (dim == 2)
	[elem2edge, isbdedge, bdnode, edge] = ...
		auxstructure(mesh, 'elem2edge','isbdedge','bdnode', 'edge');	
    NE = numel(isbdedge);
	femspace.elem2dof = [N + double(elem2edge), ...
        N + NE + double(elem2edge)];
    aux = find(isbdedge);
	femspace.bd_dof   = [bdnode; N + aux; N + NE + aux];
	femspace.ndof = N + 2*NE + NT;
    femspace.lastbdnode = N + 2*NE;
else
	[elem2edge, isbdedge, bdnode, ...
        elem2face, isbdface, edge] = ...
		auxstructure(mesh, 'elem2edge','isbdedge','bdnode', ...
        'elem2face', 'isbdface', 'edge');	
    NE = numel(isbdedge);
    NF = numel(isbdface);
	femspace.elem2dof = [N + double(elem2edge), ...
        N + NE + double(elem2edge), ...
        N + 2*NE + double(elem2face)];
    aux = find(isbdedge);
	femspace.bd_dof   = [bdnode; N + aux; N + NE + aux; ...
        N + 2*NE + find(isbdface)];
	femspace.ndof = N + 2*NE + NF;
end

if dim>1
	LNE = (dim+1)*dim/2;
	tr = get_subsimplices(1:dim+1, 1);
	for i=1:LNE
		idx = edge(elem2edge(:,i),1) ~= mesh.elem(:,tr(i,1));
        i1 = i;
        i2 = LNE + i;
        femspace.elem2dof(idx,[i1,i2]) = femspace.elem2dof(idx,[i2,i1]);
	end
end

femspace.deg = 3;
femspace.elemtype = 'p3';
if dim == 1
	femspace.nodelambda = [eye(dim+1); [2/3,1/3]; [1/3,2/3]];
elseif dim ==2
	femspace.nodelambda = [eye(dim+1); [2/3,1/3,0]; [2/3,0,1/3]; [0,2/3,1/3];
		[1/3,2/3,0]; [1/3,0,2/3]; [0,1/3,2/3]; [1/3,1/3,1/3]];
elseif (dim == 3)
	femspace.nodelambda = [eye(dim+1); 1/3*[... 
		[2,1,0,0]; [2,0,1,0]; [2,0,0,1]; [0,2,1,0]; [0,2,0,1]; [0,0,2,1]; ...
		[1,2,0,0]; [1,0,2,0]; [1,0,0,2]; [0,1,2,0]; [0,1,0,2]; [0,0,1,2]; ...
		[1,1,1,0]; [1,1,0,1]; [1,0,1,1]; [0,1,1,1];]];
else
	error('nodelambda not defined');
end

end