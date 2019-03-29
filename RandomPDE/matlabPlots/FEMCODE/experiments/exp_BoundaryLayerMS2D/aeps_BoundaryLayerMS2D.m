function out = aeps_BoundaryLayerMS2D(xmac, xmic, k, l)
%AEPS_POISSONMS_2D Summary of this function goes here
%   Detailed explanation goes here


NP = size(xmic,1);

scal = 1e-3;

% periodic
if k == l
	out = scal*(sin(2*pi*xmic(:,1))+1.125).*(cos(2*pi*xmic(:,2))+1.125)/(9*sqrt(17)/64);
else
	out = zeros(NP,1);
end

end

