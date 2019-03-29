function L2norm = get_L2error_exact( mesh, femspace, ff, f)
%GET_L2ERROR gets L2 error of two discrete functions

dim = size(mesh.elem,2) - 1;
deg = 2*femspace.deg;

[lambda, weight] = quadpts(dim, deg);
if ~isfield(mesh,'volume'), mesh.volume = simplex_volume(mesh); end

res = 0;
for i=1:numel(weight)
	res = res + weight(i) * ( ...
		evalf(mesh, 'all', femspace, lambda(i,:), ff) - ...
		evalf(mesh, 'all', [], lambda(i,:), f)).^2;
end
res = sum(bsxfun(@times, sum(res,2), mesh.volume), 1);

L2norm =sqrt(res);
end

