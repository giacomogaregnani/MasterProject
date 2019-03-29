function xx = ex3_dg(x,d)
%EXACT SOLUTION FOR DG DIFFUSION TEST IN 3D

if nargin < 2
	d=[0,0,0];
elseif d > 0
	a=[0,0,0];
	a(d) = 1;
	d=a;
end

if (sum(d)==0)
    xx =sin(pi*x(:,1)) .* sin(pi*x(:,2)) .* sin(pi*x(:,3));
elseif (sum(d) == 1) && (d(1) == 1)
    xx = pi * cos(pi*x(:,1)) .* sin(pi*x(:,2)) .* sin(pi*x(:,3));
elseif (sum(d) == 1) && (d(2) == 1)
    xx = pi * sin(pi*x(:,1)) .* cos(pi*x(:,2)) .* sin(pi*x(:,3));
elseif (sum(d) == 1) && (d(3) == 1)
    xx = pi * sin(pi*x(:,1)) .* sin(pi*x(:,2)) .* cos(pi*x(:,3));
else
	error('unimplemented derivative');
end

end

