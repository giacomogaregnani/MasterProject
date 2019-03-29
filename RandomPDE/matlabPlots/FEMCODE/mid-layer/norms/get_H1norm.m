function [ H1norm, L2norm, H1seminorm] = get_H1norm(mesh, femspace, sol)
%GET_H1NORM Summary of this function goes here
%   Detailed explanation goes here

dim = size(mesh.elem,2) - 1;
deg = 2*(femspace.deg-1);

[lambda, weight] = quadpts(dim, deg);
if ~isfield(mesh,'volume'), mesh.volume = simplex_volume(mesh); end

semi = 0;
for i=1:numel(weight)
	for j=1:dim
		semi = semi + weight(i) * ...
			evalf(mesh, 'all', femspace, lambda(i,:), sol, j).^2;
	end	
end

semi = sum(bsxfun(@times, sum(semi,2), mesh.volume), 1);

H1seminorm = sqrt(semi);
L2norm = get_L2norm(mesh, femspace, sol);
H1norm = sqrt(semi + L2norm.^2);
end

