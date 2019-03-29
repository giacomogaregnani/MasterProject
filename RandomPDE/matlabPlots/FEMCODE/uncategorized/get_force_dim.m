function [FD, SF] = get_force_dim(mesh, f)
%GET_FORCE_DIM Summary of this function goes here
%   Detailed explanation goes here

if isnumeric(f)
	s = size(f);
else
	s = size(f(mesh.node(1,:)));
end
if numel(s) <= 2
	FD = s(2); SF = 1;
elseif numel(s) >= 3
	FD = s(2); SF = prod(s(3:end));
end

