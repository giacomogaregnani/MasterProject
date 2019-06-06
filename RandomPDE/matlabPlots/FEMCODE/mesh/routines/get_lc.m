function loc_lambda = get_lc(cmesh, fmesh, father, lambda)
%GET_LC maps barycentric coordinates of a fine mesh to a coarser mesh
%
% loc_lambda = get_lc(cmesh, fmesh, father, lambda) computes the
% barycentric coordinates in a coarse mesh cmesh in respective elements to
% points given in a fine mesh fmesh by barycentric coordinates fmesh. To
% faciliate the computation, the field father describes the inclusive
% relation between the meshes.
%
%   loc_lambda(i,:) are the barycentric coordinates in the element
%   father(i) in the mesh cmesh describing the point in the element number
%   i in the mesh fmesh with barycentric coordinates lambda(i,:). If lambda
%   is just a vector, we use the same barycentric coordinates in the fine
%   mesh in every occasion.
% 
%   father(i) = j   means that the element number i in fmesh is
%   geometrically included in the element number j in cmesh.

if ~isfield(fmesh,'periodic'), fmesh.periodic = false; end 

dim = size(fmesh.node, 2);
dim_el = size(fmesh.elem, 2) - 1;
NP = size(fmesh.elem, 1);
if dim~=dim_el
    l = simplex_diameter(cmesh);
    b = get_rc(fmesh, [], 'all', [], lambda) - ...
        cmesh.node(cmesh.elem(father,1),:);
    c = sqrt(sum(b.^2,2));
    loc_lambda = [1 - c./l(father), c./l(father)];
    return
end
A = zeros(NP, dim, dim);

for i=1:dim
	A(:,:,i) = reshape(get_vec(cmesh, ...
		[cmesh.elem(father,1), cmesh.elem(father,i+1)]), [NP,dim,1]);
end

if ~fmesh.periodic
	rhs = get_rc(fmesh, [], 'all', [], lambda) - ...
		cmesh.node(cmesh.elem(father,1),:);
else
	rhs = per_diff(get_rc(fmesh, [], 'all', [], lambda) - ...
		cmesh.node(cmesh.elem(father,1),:), cmesh.box);
end
loc_lambda = solve_all(A, rhs);
loc_lambda = [1-sum(loc_lambda,2), loc_lambda];
end
