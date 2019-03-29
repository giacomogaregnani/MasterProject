function lambda2 = get_lc2(mesh, whe, x)
%GET_LC maps cartesian to barycentric coordinates of specified elements


dim = size(mesh.node, 2);
NP = numel(whe);
A = zeros(NP, dim, dim);

for i=1:dim
	A(:,:,i) = reshape(get_vec(mesh, ...
		[mesh.elem(whe,1), mesh.elem(whe,i+1)]), [NP,dim,1]);
end

if ~isfield(mesh,'periodic') || ~mesh.periodic
	rhs = x - mesh.node(mesh.elem(whe,1),:);
else
	rhs = per_diff(x - mesh.node(mesh.elem(whe,1),:), mesh.box);
end
lambda2 = solve_all(A, rhs);
lambda2 = [1-sum(lambda2,2), lambda2];
end

