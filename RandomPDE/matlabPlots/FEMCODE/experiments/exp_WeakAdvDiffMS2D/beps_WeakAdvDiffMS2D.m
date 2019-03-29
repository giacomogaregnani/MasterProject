function out = beps_WeakAdvDiffMS2D(xmac, xmic, k)
%BEPS_WEAKADVDIFFMS2D Summary of this function goes here
%   Detailed explanation goes here


% NP = size(xmic,1);

% locally periodic
if k == 1
	out = 1./(xmac(:,1).^3+3+2*sqrt(17)/8*1./(sin(2*pi*xmic(:,1))+9/8));
elseif k==2
    out = 1./(xmac(:,2).^2+0.05+(xmac(:,1).*xmac(:,2)+1)...
        .*(2*sqrt(17)/8*1./(cos(2*pi*xmic(:,2))+9/8)));
else
    error('invalid k');
end
    

end

