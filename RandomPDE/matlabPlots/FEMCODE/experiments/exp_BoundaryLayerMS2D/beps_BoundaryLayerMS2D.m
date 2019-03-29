function out = beps_BoundaryLayerMS2D(xmac, xmic, k)
%BEPS_WEAKADVDIFFMS2D Summary of this function goes here
%   Detailed explanation goes here


% NP = size(xmic,1);

% periodic
out = (sin(2*pi*xmic(:,1))+1.125).*(cos(2*pi*xmic(:,2))+1.125)/(9*sqrt(17)/64);

end

