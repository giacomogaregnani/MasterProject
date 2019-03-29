function [mesh, updatenode] = deperiodize(mesh)
%deperiodize makes non-periodic mesh from a periodic one 
%
% [mesh, updatesol] = deperiodize(mesh) takes a periodic mesh and creates a
% non-periodic equivalent, where suitable nodes are doubled (quadrupled,
% ...). The i-th element of the vector updatenode is the index of the node
% in the input mesh corresponding to the i-th node in the output mesh.

N = size(mesh.node, 1);
d = size(mesh.node, 2);
NT = size(mesh.elem, 1);
dim = size(mesh.elem, 2) -1;

if ~isfield(mesh,'periodic') || ~mesh.periodic, 
    updatenode = 1:N;
    return; 
end

tol = 1e-12;

updatenode = 1:N;
mesh.bary=get_rc(mesh);
C = zeros(NT*(dim+1),d);

for i=1:dim+1
	C((i-1)*NT+1:i*NT,:) = ...
		get_rc(mesh,[],'all',i,[],struct('adjust_to','bary'));
end

for i=1:d
	A = mesh.node(mesh.elem(:),i);
	B = C(:,i);
	nch = abs(A-B)>tol;
	[inch,~,b] = unique(mesh.elem(nch));

	N = size(mesh.node,1);
	newN = numel(inch);
	mesh.node(N+1:N+newN,[1:i-1, i+1:d]) = ...
		mesh.node(inch,[1:i-1, i+1:d]);
	mesh.node(N+1:N+newN,i) = ...
		mesh.box(2*i-1) + (mesh.box(2*i) - mesh.node(inch,i));
	mesh.elem(nch) = N + b;
	updatenode(N+1:N+newN) = updatenode(inch);
end
mesh.periodic = false;
end

