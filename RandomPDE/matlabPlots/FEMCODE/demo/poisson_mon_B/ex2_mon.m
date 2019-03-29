function xx = ex2_mon(x,d)
%EXACT SOLUTION FOR MONOTONE PROBLEM IN 2D

if nargin<2 || (sum(d)==0)
    xx = 8*sin(pi*x(:,1)).*x(:,2).*(1-x(:,2));
elseif sum(d) == 1 && d(1) == 1
    xx = 8*cos(pi*x(:,1))*pi.*x(:,2).*(1-x(:,2));
else
    xx = 8*sin(pi*x(:,1)).*(1-x(:,2))-8*sin(pi*x(:,1)).*x(:,2);
end

end

