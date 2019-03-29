function par = x2par_X(x)
par = [ cos(2*pi*(x(:,2)-x(:,1))/4)/5,   cos(2*pi*(x(:,2)+x(:,1))/4)/5];
end