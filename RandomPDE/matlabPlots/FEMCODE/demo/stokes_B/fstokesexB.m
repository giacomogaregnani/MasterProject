function F = fstokesexB( x, der )
%FSTOKES Summary of this function goes here
%   Detailed explanation goes here
%N=size(x,1);
%d=size(x,2);


if (nargin == 1) || (0 == sum(der))
%	F = zeros(size(x));

 y=x(:,2);
 x=x(:,1);
 F = [...
 		(8*x.^2.*y+8*y.^3).*cos(x.^2-y.^2) ...
 		+ 2*x.*cos(x.^2+y.^2) ...
 		- 8*y.*sin(x.^2-y.^2),...
 		(8*x.^3+8*x.*y.^2).*cos(x.^2-y.^2) ...
 		+ 2*y.*cos(x.^2+y.^2) ...
 		+ 8*x.*sin(x.^2-y.^2),...
 		];
else
	error('not defined');
end

