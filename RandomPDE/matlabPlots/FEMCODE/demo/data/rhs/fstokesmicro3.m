function out = fstokesmicro3( x, der )
%F2 Summary of this function goes here
%   Detailed explanation goes here
N=size(x,1);
d=size(x,2);
if (nargin == 1) || (0 == sum(der))
	out = [zeros(N, 2), ones(N, 1), zeros(N, d - 3)];
else
	out = zeros(N,d);
end
end
