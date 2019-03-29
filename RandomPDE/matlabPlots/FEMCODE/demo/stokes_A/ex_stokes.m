function xx = ex_stokes(x, d)
%EX_STOKES Summary of this function goes here
%   Detailed explanation goes here
y=x(:,2);
x=x(:,1);

if nargin<2 || (sum(d)==0)
    xx = [-2*(-1+x).^2.*x.^2.*(-1+y).*y.*(-1+2*y), ...
		   2*(-1+y).^2.*y.^2.*(-1+x).*x.*(-1+2*x)];
elseif sum(d) == 1 && d(1) == 1
    xx = [4*(1-2*x).*(-1+x).*x.*(-1+y).*y.*(-1+2*y), ...
		2*(1-6*x+6.*x.^2).*(-1+y).^2.*y.^2];
else
	xx = [-2*(-1+x).^2.*x.^2.*(1-6*y+6*y.^2), ...
		4*(-1+x).*x.*(-1+2*x).*(-1+y).*y.*(-1+2*y)];
end