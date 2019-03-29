function res = evalf(mesh, whe, femspace, lambda, f, der, varargin)
%EVALF evaluates a function at a set of points on a mesh

dim = size(mesh.elem, 2) - 1;

if (nargin < 6) || isempty(der)
	der = zeros(1,dim);
end

if dim>1 && numel(der) == 1
	der = der == 1:dim;
end

if isnumeric(f)
	if sum(der)>0 && ~isfield(mesh, 'dlambda'),
		mesh = gradbasis(mesh);
	end
	
	res = 0;
	for i = 1 : femspace.ldof
		res = res + ...
			bsxfun(@times,...
			f(get_dof(mesh, whe, femspace, i),:,:), ...
			evalb(mesh, whe, lambda, i, der, femspace.elemtype));
	end
else
	pxy = get_rc(mesh, [], whe, [], lambda);
	if sum(der) > 0
		res = f(pxy, der, varargin{:});
	else
		res = f(pxy, varargin{:});
	end
end
end

