function out = aeps_PoissonMS2D(xmac, xmic, k, l)
%AEPS_POISSONMS_2D Summary of this function goes here
%   Detailed explanation goes here


NP = size(xmic,1);

% periodic
% if k == l
%     out = xmac(:,1) + (sin(2*pi*xmic(:,1)).^2 + sin(2*pi*xmic(:,2)).^2);
% else
%     out = zeros(NP,1);    
% end

% locally periodic
% AFFINE
% if k == 1 && l == 1
% 	out = (16*(xmac(:,1).^2-xmac(:,1)).*(xmac(:,2).^2-xmac(:,2))+1).*(cos(2*pi*xmic(:,1)).^2 + 1);
% elseif k==2 && l == 2
%     out = (16*(xmac(:,1).^2-xmac(:,1)).*(xmac(:,2).^2-xmac(:,2))+1).*(sin(2*pi*xmic(:,2)) + 2);
% else
% 	out = zeros(NP,1);
% end

% NONAFFINE
f = @(A, x, x0, y0, sx, sy) A*exp( - ( (x(:,1) - x0).^2/(2*sx^2) + (x(:,2) - y0).^2/(2*sy^2) ));
xx = @(x) f(1, x, 2/3, 2/3, 0.15, 0.15) + 2;
if k == 1 && l == 1
    out = 300.*(((xmic(:,1)-.5).^2 + (xmic(:,2)-.5).^2)>=.5^2) + 0.1 + 0.*zeros(NP,1);%sqrt((xx(xmac) + sin(4*pi*xmic(:,1))).*(xmac(:,1).^2 + 1.2 + sin(2*pi*xmic(:,1))));
elseif k == 2 && l == 2
    out = 300.*(((xmic(:,1)-.5).^2 + (xmic(:,2)-.5).^2)>=.5^2) + 0.1 + 0.*zeros(NP,1);%sqrt((xx(xmac) + cos(5*pi*xmic(:,2))).*(xmac(:,2).^2.*cos(2*pi*xmic(:,2)) + xmac(:,1) + 1.5));
else
    out = zeros(NP,1);
end

% experiment flux
% if k == 1 && l == 1
%     out = 1./sqrt((xmac(:,1).^2 + sin(2*pi*xmic(:,1)) + 1.2).*(xmac(:,1).*xmac(:,2) + sin(4*pi*xmic(:,1)) + 1.5));
% elseif k == 2 && l == 2
%     out = 1./((xmac(:,1).*xmac(:,2) + sin(5*pi*xmic(:,2)) + 1.2).*(xmac(:,2).^2.*cos(2*pi*xmic(:,2)) + xmac(:,1) + 1.5));
% else
%     out = zeros(NP,1);
% end

end

