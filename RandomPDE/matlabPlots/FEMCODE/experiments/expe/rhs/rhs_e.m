function out = rhs_expe(p, der)
%FSTOKES_MS Summary of this function goes here
%   Detailed explanation goes here
if (nargin > 1) && (sum(der) > 0)
	error('Derivatives not implemented');
else
	d = size(p,2);
	NP = size(p,1);
	switch d
		case 1
			out = -[zeros(NP,1),  ones(NP,1)];
		case 2
			out = -[zeros(NP,2), ones(NP,1)];
		case 3
			out = -[zeros(NP,3), ones(NP,1)];
		otherwise
			error('Dimension > 3 not implemented');
	end
end

end

