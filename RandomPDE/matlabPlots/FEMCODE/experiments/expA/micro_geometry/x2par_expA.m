function par = x2par_expA(x)
par = [ cos(2*pi*(x(:,2)-x(:,1))/3)/5,   cos(2*pi*(x(:,2)+x(:,1))/3)/5];
end