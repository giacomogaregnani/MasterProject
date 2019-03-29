function mesh = cleanBdElem(mesh)
%CLEANBDELEM REMOVES ELEMENTS WITH MORE THAN ONE BOUNDARY SIDES
dim = size(mesh.node,2);

NT = size(mesh.elem);
NTold = size(mesh.elem)+1;
while sum(NT<NTold)
	isbdnode = auxstructure(mesh,'isbdnode');
	idx = (sum(isbdnode(mesh.elem),2) < dim+1);
	mesh.elem = mesh.elem(idx,:);
	NTold = NT;
	NT = size(mesh.elem);
end
end

