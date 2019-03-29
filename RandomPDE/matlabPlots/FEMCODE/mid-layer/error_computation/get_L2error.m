function [L2norm, elemL2norm] = ...
        get_L2error(cmesh, cfemspace, cf, fmesh, ffemspace, ff, father)
%GET_L2ERROR gets L2 error of two discrete functions

dim = size(cmesh.elem,2) - 1;
deg = max([2*cfemspace.deg, 2*ffemspace.deg]);

[lambda, weight] = quadpts(dim, deg);
if ~isfield(fmesh,'volume'), fmesh.volume = simplex_volume(fmesh); end

res = 0;
for i=1:numel(weight)
        loc_lambda_i = get_lc(cmesh, fmesh, father, lambda(i,:));
        res = res + weight(i) * ( ...
                evalf(fmesh, 'all', ffemspace, lambda(i,:), ff) - ...
                evalf(cmesh, father, cfemspace, loc_lambda_i, cf)).^2;
end
res = bsxfun(@times, sum(res,2), fmesh.volume);

elemL2norm = sqrt(res);
L2norm = sqrt(sum(res, 1));
end
