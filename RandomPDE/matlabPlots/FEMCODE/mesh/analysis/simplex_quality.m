function q = simplex_quality(mesh)
%SIMPLEX_QUALITY computes the quality of elements in a mesh
% 
% q = simplex_quality(mesh) computes the quality index based on the ratio
% of the inscribed and the subscribed circle (sphere) of all elements in a
% mesh. The quality index is 1 for equilateral triangle or for a regular
% tetrahedron. It decreases to 0 with singular deformations.

dim = size(mesh.elem,2) - 1;
NT = size(mesh.elem,1);
d = size(mesh.node,2);

if (dim == 1) 
	q=ones(1,size(mesh.elem,1)); 
elseif (dim == 2) 
	% well-known formula
	a = sqrt(sum(get_vec(mesh,[],'all',[1,2]).^2,2)); % sizes of edges
	b = sqrt(sum(get_vec(mesh,[],'all',[1,3]).^2,2));
	c = sqrt(sum(get_vec(mesh,[],'all',[2,3]).^2,2));
	S = simplex_volume(mesh); % area
	r = 2*S ./(a+b+c);        % inradius
	R = a.*b.*c ./ (4*S);     % circumradius
	q = 2*r ./ R;             % quality
elseif (dim == 3)
	% not so well-known formula
	d12 = get_vec(mesh,[],'all',[1,2]); % edges as vectors
	d13 = get_vec(mesh,[],'all',[1,3]);
	d14 = get_vec(mesh,[],'all',[1,4]);
	d23 = get_vec(mesh,[],'all',[2,3]);
	d24 = get_vec(mesh,[],'all',[2,4]);
	d34 = get_vec(mesh,[],'all',[3,4]);
	V   = simplex_volume(mesh); % volume
	s1  = sqrt(sum(cross(d12,d13,2).^2,2))/2; % areas of faces
	s2  = sqrt(sum(cross(d12,d14,2).^2,2))/2;
	s3  = sqrt(sum(cross(d13,d14,2).^2,2))/2;
	s4  = sqrt(sum(cross(d23,d24,2).^2,2))/2;
	p1  = sqrt(sum(d12.^2,2)).*sqrt(sum(d34.^2,2)); % 
	p2  = sqrt(sum(d23.^2,2)).*sqrt(sum(d14.^2,2));
	p3  = sqrt(sum(d13.^2,2)).*sqrt(sum(d24.^2,2));
	q   = 216*V.^2./(s1+s2+s3+s4) ./ sqrt((p1+p2+p3) .* ...
		(p1+p2-p3) .* (p1+p3-p2) .* (p2+p3-p1)); % quality
elseif (d==dim) 
%% GENERAL ALGORITHM for dim == d to compute the circumradius
% Let:		A1, A2, ..., Adim+1 be vertices of a simplex and
% denote:	Vi = (A(i+1) - A1) / 2
% Denote the circumcenter by O and let us search V := O - A1 
% For every i = 1:d, the point O must lie at the same distance to the
% points A1 and A(i+1) and so it must lie on a hyperplane orthogonal to Vi
% passing through the midpoint of the edge A1 -- A(i+1). This can be simply
% expressed as
%					Vi dot V = (Vi dot Vi)
% These relations are necessary and sufficient to determine a unique
% solution V and R = |V|.
	A=zeros(NT, d, d);
	rhs = zeros(NT, d);
	for i=1:dim
		vi = get_vec(mesh, [], 'all', [1,i+1])/2;
		A(:,i,:) = reshape(vi, [NT, 1, dim]);
		rhs(:,i) = sum(vi.^2, 2);
	end
	sol = solve_all(A,rhs);
	R = sqrt(sum(sol.^2,2));				% CIRCUMRADIUS
	V = simplex_volume(mesh);				% VOLUME
	subsim = subsimplex_volume(mesh, 1);	% LENGTHS OF SUBSIMPLICES		
	r = dim * V ./ sum(subsim, 2);			% INRADIUS
	q = dim * r ./ R;						% QUALITY 
else 
	error('Case (d > dim > 3) not implemented.');
end
end