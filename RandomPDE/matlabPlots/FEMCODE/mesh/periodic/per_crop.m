function p = per_crop(p, box)
%PER_CROP takes points from R^d and finds their corresp. coor. on a torus
%
% p = per_crop(p, box) takes the points in R^d from the  array p (one row 
% is one point) and finds where they are mapped to the torus given by the
% array box. 

d = size(p,2);
dimensions = box(2:2:end) - box(1:2:end); % sizes of domain

for i=1:d
	ind = ~( (p(:,i) >= box(2*i-1)) & (p(:,i) < box(2*i)) );
	p(ind,i) = p(ind,i) - dimensions(i) * ...
		floor((p(ind,i) - box(2*i-1)) / dimensions(i));
end
end