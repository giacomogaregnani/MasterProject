function  res  = get_L2scalar(mesh, femspace, sol1, sol2)
%GET_L2NORM Summary of this function goes here
%   Detailed explanation goes here

dim = size(mesh.elem,2) - 1;
deg = 2*femspace.deg;

if ~isfield(mesh,'volume'), mesh.volume = simplex_volume(mesh); end
[lambda, weight] = quadpts(dim, deg);
	
res = 0; %weight(1) * evalf(mesh, [], femspace, lambda(1,:), f).^2;
for i=1:numel(weight)
	res = res + weight(i) * ...
		evalf(mesh, 'all', femspace, lambda(i,:), sol1) .* ...
		evalf(mesh, 'all', femspace, lambda(i,:), sol2);
end
res = sum(bsxfun(@times, sum(res,2), mesh.volume), 1);
end

