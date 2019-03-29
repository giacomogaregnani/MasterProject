function res = p1_b(mesh, whe, lambda, bf, der)
%P1_B Summary of this function goes here
%   Detailed explanation goes here
if strcmp(whe,'all')
	NT = size(mesh.elem,1);
else
	NT = numel(whe);
end

if (sum(der)==0)
	res = lambda(:,bf);
	if size(res,1) == 1
		res = repmat(res, [NT, 1] );
	end
elseif (sum(der)==1)
	if strcmp(whe,'all')
		res = mesh.dlambda(:, der > 0, bf);
	else
		res = mesh.dlambda(whe, der > 0, bf);
	end
else
	res = zeros(NT,1);
end

end
