function [femspace, res] = restore(mesh, f, elemtype)
%RESTORE reconstructs a discontinuous discrete function from quad. points
%
% takes a function that is 

% TODO: here we implicitly use the knowledge of quadrature formulas and of
% the degrees of freedom numbers. This should be made general.

dim = size(mesh.elem, 2) - 1;

if strcmp(elemtype,'p1') || strcmp(elemtype,'p1d')
	femspace = get_femspace(mesh,'p0d');
	res = f;
    return;
elseif strcmp(elemtype,'p2') || strcmp(elemtype,'p2d')
	femspace = get_femspace(mesh,'p1d');

	alpha = (1 - 1/sqrt(dim+2)) / (dim+1);
	beta  = 1 - dim*alpha;
	
	fsum = sum(f,4);
    s = size(f);
    
	res = bsxfun(@plus,f,-alpha*fsum) / (beta-alpha);
	res = permute(res, [4,1,2,3]);
	res = reshape(res, [], s(2), s(3));	
elseif (dim == 2) && (strcmp(elemtype,'p3') || strcmp(elemtype,'p3d'))
	femspace = get_femspace(mesh,'p2d');
	
	lambda = quadpts(2,4);
	a = lambda(1,1);
	b = lambda(1,2);
	c = lambda(4,1);
	d = lambda(4,2);
	M = [...
		2*a*(a-0.5), 2*b*(b-0.5), 2*b*(b-0.5), 4*a*b, 4*a*b, 4*b*b;
		2*b*(b-0.5), 2*a*(a-0.5), 2*b*(b-0.5), 4*a*b, 4*b*b, 4*a*b;
		2*b*(b-0.5), 2*b*(b-0.5), 2*a*(a-0.5), 4*b*b, 4*a*b, 4*a*b;
		2*c*(c-0.5), 2*d*(d-0.5), 2*d*(d-0.5), 4*c*d, 4*c*d, 4*d*d;
		2*d*(d-0.5), 2*c*(c-0.5), 2*d*(d-0.5), 4*c*d, 4*d*d, 4*c*d;
		2*d*(d-0.5), 2*d*(d-0.5), 2*c*(c-0.5), 4*d*d, 4*c*d, 4*c*d;
		];
	M = inv(M);
	s = size(f);
	res = zeros(size(f));
    
	for i=1:6
		for j=1:6
			res(:,:,:,i) = res(:,:,:,i) + ...
				M(i,j) * f(:,:,:,j);
		end
	end
	res = permute(res, [4,1,2,3]);
	res = reshape(res, [], s(2), s(3));
else
	error('not yet implemented');
end
end

