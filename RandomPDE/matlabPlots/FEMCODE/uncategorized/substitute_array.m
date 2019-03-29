function A = substitute_array(A, vfrom, vto)
%SUBSTITUTE_ARRAY makes a substitution in an array
%   Detailed explanation goes here
notupdated = true(size(A));
for i=1:numel(vfrom)
	update = notupdated & (A==vfrom(i)); 
	A(update) = vto(i);
	notupdated = notupdated & ~update;
end
end

