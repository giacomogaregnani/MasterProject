function int = get_integral(mesh, femspace, f, varargin)
%GET_AVERAGE computes an integral of a discrete or continuous function

dim = size(mesh.elem,2) - 1;
if ~isfield(mesh,'volume'), mesh.volume = simplex_volume(mesh); end

[lambda, weight] = quadpts(dim, femspace.deg);

res = 0;
for k=1:numel(weight)
	res = res + weight(k) * ...
        evalf(mesh, 'all', femspace, lambda(k,:), f, [], varargin{:});
end
int = sum(bsxfun(@times, res, mesh.volume), 1);
end

