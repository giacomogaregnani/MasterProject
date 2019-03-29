function [H1norm, L2norm, H1seminorm] = ...
	get_H1error_exact(mesh, femspace, ff, f)
%GET_H1ERROR gets H1, L2 and H1seminorm errors of two discrete functions

dim = size(mesh.elem,2) - 1;
deg = 2*femspace.deg;

[lambda, weight] = quadpts(dim, deg);
if ~isfield(mesh,'volume'), mesh.volume = simplex_volume(mesh); end

semi = 0;
for i=1:numel(weight)
	for j=1:dim
		semi = semi + weight(i) * ( ...
			evalf(mesh, 'all', femspace, lambda(i,:), ff, j) - ...
			evalf(mesh, 'all', [], lambda(i,:), f, j)).^2;
	end	
end
semi = sum(bsxfun(@times, sum(semi,2), mesh.volume), 1);

H1seminorm = sqrt(semi);
L2norm = get_L2error_exact(mesh, femspace, ff, f);
H1norm = sqrt(semi + L2norm.^2);
end

