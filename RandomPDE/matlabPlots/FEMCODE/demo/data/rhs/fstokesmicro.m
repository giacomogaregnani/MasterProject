function [ out] = fstokesmicro( x, der )
%FSTOKESMICRO Summary of this function goes here
%   Detailed explanation goes here

N=size(x,1);
d=size(x,2);
if (nargin == 1) || (0 == sum(der))
	out = zeros(N,d,d);
	for i=1:d
		out(:,i,i) = 1;
	end
else
	out = zeros(N,d,d);
end

end

