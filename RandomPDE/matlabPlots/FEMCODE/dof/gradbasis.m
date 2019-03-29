function mesh = gradbasis(mesh)
%GRADBASIS returns gradients of barycentric coordinates in each element
%
% mesh = gradbasis(mesh) adds a field mesh.dlambda which is an array of
% size NT x dim x dim+1. Number mesh.dlambda(n, m, lam) is the derivative
% with respect to x_m of the lam-th barycentric coordinate in the element
% number m.
%

if isfield(mesh,'dlambda'), 
	if ~isfield(mesh,'volume')
		mesh.volume = simplex_volume(mesh);
	end
	return;
end

dim = size(mesh.elem,2)-1;
fac = factorial(dim);
NT = size(mesh.elem,1);

oriented_volume = simplex_volume(mesh, true);
if (dim == 1)
	mesh.dlambda(1:NT,:,1) = -1./oriented_volume;
	mesh.dlambda(1:NT,:,2) = 1./oriented_volume;	
elseif (dim == 2)
	ve1 = get_vec(mesh,[],'all',[2,3]);
	ve2 = get_vec(mesh,[],'all',[3,1]);
	ve3 = get_vec(mesh,[],'all',[1,2]);
	if ~isfield(mesh,'dlambda')
		mesh.dlambda(1:NT,:,3) = ...
			[-ve3(:,2), ve3(:,1)] ./ repmat(fac*oriented_volume,[1, dim]);
		mesh.dlambda(1:NT,:,1) = ...
			[-ve1(:,2), ve1(:,1)] ./ repmat(fac*oriented_volume,[1, dim]);
		mesh.dlambda(1:NT,:,2) = ...
			[-ve2(:,2), ve2(:,1)] ./ repmat(fac*oriented_volume,[1, dim]);
	end
	
elseif (dim == 3)
		%% SPECIAL ORIENTATION
	face = uint32([mesh.elem(:,[2 4 3]); mesh.elem(:,[1 3 4]); ...
		mesh.elem(:,[1 4 2]); mesh.elem(:,[1 2 3])]);
	v12 = get_vec(mesh,face(:,[1 2])); 
	v13 = get_vec(mesh,face(:,[1 3])); 
	normal = cross(v12,v13,2);
	if ~isfield(mesh,'dlambda')
		mesh.dlambda(1:NT,:,4) = ...
			normal(3*NT+1:4*NT,:) ./ repmat(fac*oriented_volume,[1,dim]);
		mesh.dlambda(1:NT,:,1) = ...
			normal(1:NT,:) ./        repmat(fac*oriented_volume,[1,dim]);
		mesh.dlambda(1:NT,:,2) = ...
			normal(NT+1:2*NT,:) ./   repmat(fac*oriented_volume,[1,dim]);
		mesh.dlambda(1:NT,:,3) = ...
			normal(2*NT+1:3*NT,:) ./ repmat(fac*oriented_volume,[1,dim]);
	end
else
	error('General dimension gradbasis not yet implemented');
end
if ~isfield(mesh,'volume')
	mesh.volume = abs(oriented_volume);
end
end

