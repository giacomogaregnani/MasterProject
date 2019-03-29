function vec = per_diff(vec, box)
%PER_DIFF normalizes vectors given on a torus to their minimal size
%
% vec = per_diff(vec, box) takes the array vec containing vectors 
% (one row = one vector) that are small on a torus given by the array box. 
% This function finds equivalent vectors with minimal size.

d = size(vec,2);
dimensions = box(2:2:end) - box(1:2:end); % sizes of domain

for i=1:d
	ind = ~(abs(vec(:,i)) < dimensions(i)/2); % find where the i'th 
	vec(ind,i) = vec(ind,i) - dimensions(i) * ...
		round(vec(ind,i) / dimensions(i));
end
end