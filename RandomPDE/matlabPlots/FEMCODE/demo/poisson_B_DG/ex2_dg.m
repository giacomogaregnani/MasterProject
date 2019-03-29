function xx = ex2_dg(x,d)
%EXACT SOLUTION FOR DG DIFFUSION TEST IN 2D

if nargin<2 || (sum(d)==0)
    xx =sin(pi*x(:,1)) .* sin(pi*x(:,2));
elseif sum(d) == 1 && d(1) == 1
    xx = pi * cos(pi*x(:,1)) .* sin(pi*x(:,2));
else
    xx = pi * sin(pi*x(:,1)) .* cos(pi*x(:,2));
end

end

