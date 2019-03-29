function xx = ex_stokesB(x, d)
%EX_STOKES Summary of this function goes here
%   Detailed explanation goes here
y=x(:,2);
x=x(:,1);

if nargin<2 || (sum(d)==0)
%	xx = ones(size(x,1),2);

%xx = [y.^2-1, x.^2-1];

     xx = [2*y.*cos(x.^2-y.^2), 2*x.*cos(x.^2-y.^2)];
elseif sum(d) == 1 && d(1) == 1
%	xx = zeros(size(x,1),2);	

%xx = [zeros(size(x)), 2*x];

     xx = [-4*x.*y.*sin(x.^2-y.^2), ...
 		2*cos(x.^2-y.^2) - 4*x.^2.*sin(x.^2-y.^2)];
else
%	xx = zeros(size(x,1),2);
	
%	xx = [2*y,zeros(size(x))];

 	xx = [2*cos(x.^2-y.^2) + 4*y.^2.*sin(x.^2-y.^2), ...
 		4*x.*y.*sin(x.^2-y.^2)];
end