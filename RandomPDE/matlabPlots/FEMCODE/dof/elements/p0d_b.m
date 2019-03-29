function res = p0d_b(mesh, whe, ~, ~, der)
%P0D_B Summary of this function goes here
%   Detailed explanation goes here
if strcmp(whe,'all')
	NT = size(mesh.elem,1);
else
	NT = numel(whe);
end

if (sum(der)==0)
	res = ones(NT, 1);
else
	res = zeros(NT,1);
end
end

