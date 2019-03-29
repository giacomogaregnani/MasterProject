function out = fstokes3( x, der )
%F2 Summary of this function goes here
%   Detailed explanation goes here
N=size(x,1);
d=size(x,2);
if (nargin == 1) || (0 == sum(der))
	out = [zeros(N, 1), ones(N, 1), zeros(N, d - 2)];
else
	out = zeros(N,d);
end
end

