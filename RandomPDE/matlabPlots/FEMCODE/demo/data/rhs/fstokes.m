function [ out ] = fstokes( x, der )
%FSTOKES Summary of this function goes here
%   Detailed explanation goes here
N=size(x,1);
d=size(x,2);
if (nargin == 1) || (0 == sum(der))
	out = -[ones(N, 1), zeros(N, d - 1)];
else
	out = zeros(N,d);
end

end

