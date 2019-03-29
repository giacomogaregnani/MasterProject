function  out = repval(mat, factors)
%REPVAL substitutes each entry of Nd array by a const. matrix of size dims 
% 
% out = repval(mat, dims) takes Nd array mat and substitutes each its entry
% by a constant matrix of size factors with the same value. 

ND = numel(factors);

if ND < numel(size(mat))
	error('Repeating factor not specified for every dimension');
end

matsize=size(mat); matsize(end+1:ND) = 1;
final = matsize.*factors;

out = mat;
for i=1:ND % repval in i-th dimension
	out = reshape(out, prod(final(1:i-1)), []);
	out = repmat(out, [factors(i), ones(1,ND-1)]);
end
out = reshape(out, final);
end

%% EXAMPLE
%
% repval([1,2,3; 4,5,6],[4,5])
% 
% 1     1     1     1     1     2     2     2     2     2     3     3     3     3     3
% 1     1     1     1     1     2     2     2     2     2     3     3     3     3     3
% 1     1     1     1     1     2     2     2     2     2     3     3     3     3     3
% 1     1     1     1     1     2     2     2     2     2     3     3     3     3     3
% 4     4     4     4     4     5     5     5     5     5     6     6     6     6     6
% 4     4     4     4     4     5     5     5     5     5     6     6     6     6     6
% 4     4     4     4     4     5     5     5     5     5     6     6     6     6     6
% 4     4     4     4     4     5     5     5     5     5     6     6     6     6     6
