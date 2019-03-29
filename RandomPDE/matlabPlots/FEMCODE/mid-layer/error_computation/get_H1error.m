function [ H1norm, L2norm, H1seminorm, ...
	eH1norm, eL2norm, eH1seminorm] = ...
	get_H1error( cmesh, cfemspace, cf, fmesh, ffemspace, ff, father )
%GET_H1ERROR gets H1, L2 and H1seminorm errors of two discrete functions

dim = size(cmesh.elem,2) - 1;
deg = max([2*(cfemspace.deg-1), 2*(ffemspace.deg-1)]);

[lambda, weight] = quadpts(dim, deg);
if ~isfield(fmesh,'volume'), fmesh.volume = simplex_volume(fmesh); end

semi = 0;
for i=1:numel(weight)
	loc_lambda_i = get_lc(cmesh, fmesh, father, lambda(i,:));
	for j=1:dim
		semi = semi + weight(i) * ( ...
			evalf(fmesh, 'all', ffemspace, lambda(i,:), ff, j) - ...
			evalf(cmesh, father, cfemspace, loc_lambda_i, cf, j)).^2;
	end	
end
semi = bsxfun(@times, sum(semi,2), fmesh.volume);

eH1seminorm = sqrt(semi);
H1seminorm = sqrt(sum(semi, 1));

[L2norm, eL2norm] = ...
	get_L2error(cmesh, cfemspace, cf, fmesh, ffemspace, ff, father);

H1norm = sqrt(sum(semi, 1) + L2norm.^2);
eH1norm = sqrt(semi + eL2norm.^2);
end

