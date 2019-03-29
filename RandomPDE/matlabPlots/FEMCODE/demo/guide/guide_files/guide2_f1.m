function rhs = guide2_f1(x)
rhs = [2*(1-2*x(:,1)).^2, -x(:,2), x(:,1)];
end