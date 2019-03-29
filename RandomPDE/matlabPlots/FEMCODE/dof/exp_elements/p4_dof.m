function [femspace, mesh] = p4_dof(mesh)
%% DOF_P2 dof structure for P2 element.

dim = size(mesh.elem,2) - 1;
N  = size(mesh.node,1);
NT = size(mesh.elem,1);

femspace.ldof = ...
    (dim+1) + ...
    3 * (dim+1)*dim/2 + ...
    3 * (dim+1)*dim*(dim-1)/6 + ...
    1 * (dim+1)*dim*(dim-1)*(dim-2)/24;

if (dim == 1)
	bdnode = auxstructure(mesh,'bdnode');
	femspace.elem2dof = [];
	femspace.bd_dof = bdnode;
	femspace.ndof = N + 3*NT;
    femspace.lastbdnode = N;
elseif (dim == 2)
	[elem2edge, isbdedge, bdnode, edge] = ...
		auxstructure(mesh, 'elem2edge','isbdedge','bdnode');	
    NE = numel(isbdedge);
	femspace.elem2dof = [N + double(elem2edge), ...
        N + NE + double(elem2edge), ...
        N + 2*NE + double(elem2edge)];
    aux = find(isbdedge);
	femspace.bd_dof   = [bdnode; N + aux; N + NE + aux; N + 2 * NE + aux];
	femspace.ndof = N + 3*NE + 3*NT;
    femspace.lastbdnode = N + 3*NE;
else
	[elem2edge, isbdedge, edge, bdnode, ...
        elem2face, isbdface, face] = ...
		auxstructure(mesh, 'elem2edge','isbdedge', 'edge', 'bdnode', ...
        'elem2face', 'isbdface', 'face');	
    NE = numel(isbdedge);
    NF = numel(isbdface);
	femspace.elem2dof = [N + double(elem2edge), ...
        N + NE + double(elem2edge), ...
        N + 2*NE + double(elem2edge), ...
        N + 3*NE + double(elem2face), ...
        N + 3*NE + NF + double(elem2face), ...
        N + 3*NE + 2*NF + double(elem2face)];
    aux = find(isbdedge);
    aux2 = find(isbdface);
	femspace.bd_dof   = [bdnode; ...
        N + aux; N + NE + aux; N + 2*NE + aux; ...
        N + 3*NE + aux2; N + 3*NE + NF + aux2; N + 3*NE + 2*NF + aux2];
	femspace.ndof = N + 3*NE + 3*NF + NT;
    femspace.lastbdnode = N + 3*NE + 3*NF;
end

if dim>1
	LNE = (dim+1)*dim/2;
	mesh.s1 = ones(NT,LNE,2);
	tr = get_subsimplices(1:dim+1, 1);
	for i=1:LNE
		idx = edge(elem2edge(:,i),1) == mesh.elem(:,tr(i,1));
		mesh.s1(~idx,i,1) = 2; 
		mesh.s1(idx,i,2) = 2;
	end
end

if dim>2
	LNF = (dim+1)*dim*(dim-1)/6;
	mesh.s2 = ones(NT,LNF,3);
	tr = get_subsimplices(1:dim+1, 2);
	for i=1:LNF
		face(elem2face(:,i),1) == mesh.elem(:,tr(i,1))
	end
end

femspace.deg = 4;
femspace.elemtype = 'p4';
end