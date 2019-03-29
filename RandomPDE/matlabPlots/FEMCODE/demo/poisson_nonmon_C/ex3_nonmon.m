function xx = ex3_nonmon(x,d)
%EXACT SOLUTION FOR NONMONOTONE PROBLEM IN 3D

if nargin<2 || (sum(d)==0)
    xx = 8*sin(pi*x(:,1)).*x(:,2).*(1-x(:,2)).*sin(pi*x(:,3));
elseif sum(d) == 1 && d(1) == 1
    xx = 8*cos(pi*x(:,1))*pi.*x(:,2).*(1-x(:,2)).*sin(pi*x(:,3));
elseif (sum(d) == 1 && d(2) == 1)
    xx = (8*sin(pi*x(:,1)).*(1-x(:,2))-8*sin(pi*x(:,1)).*x(:,2)).*sin(pi*x(:,3));
elseif (sum(d) == 1 && d(3) == 1)
    xx = 8*sin(pi*x(:,1)).*x(:,2).*(1-x(:,2)).*pi.*cos(pi*x(:,3));

else
    error('wrong arguments');
end

end

