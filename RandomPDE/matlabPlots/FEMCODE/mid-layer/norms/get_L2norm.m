function  L2norm  = get_L2norm( mesh, femspace, f )
%GET_L2NORM Summary of this function goes here
%   Detailed explanation goes here

dim = size(mesh.elem,2) - 1;
deg = 2*femspace.deg;

if ~isfield(mesh,'volume'), mesh.volume = simplex_volume(mesh); end
[lambda, weight] = quadpts(dim, deg);
	
res = 0;
for i=1:numel(weight)
	res = res + weight(i) * evalf(mesh, 'all', femspace, lambda(i,:), f).^2;
end
res = sum(bsxfun(@times, sum(res,2), mesh.volume), 1);
L2norm = sqrt(res);
end

