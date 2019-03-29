function volume = simplex_volume(mesh, signed)
%SIMPLEX_VOLUME gives (signed or not) measure of elements in a mesh
%
%   v = simplex_volume(mesh) computes positive measure of simplicies in the
%   triangulation given by mesh. If the dimension of simplices is less
%   than the dimension of the space, the measure is always positive.
%
%   v = simplex_volume(mesh, signed) gives signed measure if signed = true

dim = size(mesh.elem, 2) - 1;
d   = size(mesh.node, 2);
NT  = size(mesh.elem, 1);

% ASSEMBLE EDGE VECTORS
vec = zeros(NT, d, dim);
for i=1:dim
	vec(:,:,i) = get_vec(mesh, [], 'all', [1,i+1]);
end

% VOLUME = GENERALIZED DETERMINANT / DIMENSION FACTOR
volume = detn(vec) / factorial(dim);
if ~((nargin==2) && signed)
	volume = abs(volume);
end
end
