function res = get_H1seminorm(mesh, femspace, sol)
%GET_H1NORM Summary of this function goes here
%   Detailed explanation goes here

dim = size(mesh.elem,2) - 1;

if ~isfield(mesh,'volume'), mesh.volume = simplex_volume(mesh); end
[lambda, weight] = quadpts(dim, 2*(femspace.deg-1));

res = 0;
for j=1:dim
	for i=1:numel(weight)
		res = res + weight(i) * ...
			evalf(mesh, 'all', femspace, lambda(i,:), sol, j).^2;
	end
end
res = sum(bsxfun(@times, sum(res,2), mesh.volume), 1);
res = sqrt(res);
end