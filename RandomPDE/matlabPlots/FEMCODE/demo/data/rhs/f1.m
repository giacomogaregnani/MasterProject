function [ out ] = f1( x, der )
%F1 Summary of this function goes here
%   Detailed explanation goes here

if (nargin == 1) || (0 == sum(der))
	out = ones(size(x,1), 1);
else
	out = zeros(size(x,1), 1);
end
	
end

