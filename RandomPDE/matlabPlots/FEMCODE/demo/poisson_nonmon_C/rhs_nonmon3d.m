function [ out ] = rhs_nonmon3d(x, der)
%% RHS FUNCTION FOR NONLINEAR NONMONOTONE PROBLEM

if (nargin == 1) || (0 == sum(der))
    t1 = pi .* x(:,1);
    t2 = sin(t1);
    t4 = 1 - x(:,2);
    t6 = pi .* x(:,3);
    t7 = sin(t6);
    t8 = x(:,2) .* t4 .* t7;
    t10 = 0.8e1 .* pi .* t2 .* t8;
    t11 = sin(t10);
    t12 = cos(t10);
    t14 = pi .^ 2;
    t16 = cos(t1);
    t18 = t4 .* t7;
    t35 = t2 .* x(:,2);
    t39 = t2 .^ 2;
    t40 = x(:,2) .^ 2;
    t42 = t4 .^ 2;
    t43 = t7 .^ 2;
    t52 = atan(0.8e1 .* t35 .* t18);
    t59 = cos(t6);
    t60 = t4 .* t59;
    out = -0.8e1 .* (t11 + 0.8e1 .* x(:,1) .* t12 .* t14 .* t16 .* x(:,2)...
        .* t18) .* t16 .* pi .* t8 + 0.8e1 .* (0.1e1 + x(:,1) .* t11) .* t14...
        .* t2 .* t8 - 0.64e2 .* (t2 .* t4 .* t7 - t35 .* t7) .^ 2 ./ (0.1e1...
        + 0.64e2 .* t39 .* t40 .* t42 .* t43) + 0.16e2 .* (0.2e1 + t52) .* t2...
        .* t7 - 0.8e1 .* (t12 - 0.8e1 .* x(:,3) .* t11 .* t14 .* t35 .* t60)...
        .* t2 .* x(:,2) .* t60 .* pi + 0.8e1 .* (0.1e1 + x(:,3) .* t12) .* t2...
        .* x(:,2) .* t18 .* t14;

else 
    out = zeros(size(x,1), 1);
end

end
