function out = aeps_WeakAdvDiffMS2D(xmac, xmic, k, l)
%AEPS_POISSONMS_2D Summary of this function goes here
%   Detailed explanation goes here


NP = size(xmic,1);

% locally periodic
if k == 1 && l == 1
	out = (xmac(1)^3 + 3 + 2*sqrt(17)./(8*sin(2*pi*xmic(:,1))+9)).^(-1);
elseif k==2 && l == 2
    out = (xmac(2)^2 + 0.05 + (xmac(1)*xmac(2) + 1)*2*sqrt(17)...
            ./(8*cos(2*pi*xmic(:,2))+9)).^(-1);
else
	out = zeros(NP,1);
end

end

