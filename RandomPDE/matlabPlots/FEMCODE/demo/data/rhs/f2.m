function [ out ] = f2( x, der )
%F1 Summary of this function goes here
%   Detailed explanation goes here
NP = size(x,1);
if (nargin == 1) || (0 == sum(der))
	out = ones(NP, 1, 2);
	out(:,1,2) = 20*x(:,1);
else
	out = zeros(NP, 1, 2);
end
	
end

