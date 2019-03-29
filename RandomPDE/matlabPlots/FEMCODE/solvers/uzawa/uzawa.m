function [u,p,err] = uzawa(A,B,f,g,r,omega)
%% UZAWA augmented Uzawa method for solving Stokes equations
% 
%      |A   B| |u|  = |f|
%      |B'  0| |p|  = |0|
%
% The parameter r is used in the augmented uzawa method. 
%
%      |A+rB'B  B'| |u|  = |f|
%      |B       0 | |p|  = |0|
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

Np = size(B,1);
NF = size(f,3);
if (nargin <= 3) || (g == 0)
	g = zeros(Np, 1, NF);
end
if nargin<=4
    r = 0;    omega = 2;
elseif nargin <=5;
    omega = max(r+1,2);
end

p = zeros(Np,1,NF);
u = zeros(length(f),1,NF);

A = A + r*(B'*B);

for i=1:NF
	ff = f(:,1,i) + r*(B'*g(:,1,i));
	
	err = 1;
	k = 1;
	maxIt = 100;
	while (err>1e-6) && (k<maxIt)
		uold = u(:,1,i);
		u(:,1,i) = A\(ff - B'*p(:,1,i));
		p(:,1,i) = p(:,1,i) + omega*(B * u(:,1,i) - g(:,1,i));
		err = sqrt((u(:,1,i)-uold)' * A * (u(:,1,i)-uold)) / ...
			sqrt(u(:,1,i)' * A * u(:,1,i));
		fprintf('#dof: %8.0u, Uzawa iter: %2.0u, err = %12.8g\n',...
			size(u,1)+size(p,1), k, err);
		k = k+1;
	end
end