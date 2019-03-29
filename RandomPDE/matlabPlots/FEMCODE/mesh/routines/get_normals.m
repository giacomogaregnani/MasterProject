function [normals, subsim, subsim2elem] = get_normals(mesh, whe, opver)
%GET_NORMALS computes unit normals to all subsims 
%
% normals = get_normals(mesh) generates all the subsims (in 2D - edges, in
% 3D - faces) in the lexicographical order. Then, to each of them it
% computes the outer unit normal as the subsims are belonging to the 
% element given by subsim2elem(:,3).
%
% [normals, subsim, subsim2elem, elem2subsim] = get_normals(mesh) returns
% additional arrays returned by the auxstructure function.

NT  = size(mesh.elem,1);
dim = size(mesh.elem,2) - 1;

if nargin ==1
	[subsim, subsim2elem] = auxstructure(mesh,'subsim','subsim2elem');
	
	% opposite_node is a NE x 1 vector, where opposite_node(i) is the number of
	% node that is in the element subsim2elem(i,1) but not in subsim(i,:).
	opposite_node = mesh.elem(NT*(subsim2elem(:,3)-1) + subsim2elem(:,1));
else
	subsim = mesh.elem(whe,[1:opver-1, opver+1:dim+1]);
	opposite_node = mesh.elem(NT*(opver-1) + whe);
end

NE = size(subsim,1);

% for each subsim(j,:) we solve a system - we find a vector that is 
A = zeros(NE, dim, dim);
for i=1:dim-1;
% * orthogonal to all of its basis vectors (between nodes subsim(j,1) and 
% subsim(j,i+1), i=1:dim-1)
	A(:,i,:) = reshape(get_vec(mesh,subsim(:,[1,i+1])), [NE,1,dim]);
end
rhs = zeros(NE, dim - 1);
% * and that gives negative scalar product when multiplied by the vector 
% between subsim(j,1) and opposite_node(j)
A(:,dim,:)=reshape(get_vec(mesh, [subsim(:,1), opposite_node]),[NE,1,dim]);
rhs = [rhs, repmat(-1, [NE, 1])];

% these vectors are normals and are normalized to be of unit length
normals = solve_all(A, rhs);
normals = normals ./ repmat(sqrt(sum(normals.^2,2)) ,[1,dim]); 
end

