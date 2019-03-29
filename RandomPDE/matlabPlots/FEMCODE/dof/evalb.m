function res = evalb(mesh, whe, lambda, bf, der, elemtype)
%EVALB Summary of this function goes here
%   Detailed explanation goes here

dim = size(mesh.elem,2) - 1;
h = str2func([elemtype,'_b']);

if (dim > 1) && (numel(der) == 1)
	der = der == 1:dim;
end

res = h(mesh, whe, lambda, bf, der);
end