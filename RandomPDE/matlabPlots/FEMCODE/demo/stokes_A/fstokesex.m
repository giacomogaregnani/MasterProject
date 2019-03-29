function F = fstokesex( x, der )
%FSTOKES Summary of this function goes here
%   Detailed explanation goes here
N=size(x,1);
d=size(x,2);

y=x(:,2);
x=x(:,1);

if (nargin == 1) || (0 == sum(der))
	F = -[-4*(x.^3.*(6-12*y) + x.^4.*(-3+6*y) + y.*(1-3*y+2*y.^2) - ...
		6*x.*y.*(1-3*y+2*y.^2) + 3*x.^2.*(-1+4*y-6*y.^2+4*y.^3)), ...
		4*(-3*(-1+y).^2.*y.^2 - 3*x.^2.*(1-6*y+6*y.^2) + ...
		2*x.^3.*(1-6*y+6*y.^2) + x.*(1-6*y+12*y.^2-12*y.^3+6*y.^4)) ...
		];
else
	error('not defined');

end

