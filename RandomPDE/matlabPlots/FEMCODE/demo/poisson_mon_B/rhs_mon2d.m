function [ out ] = rhs_mon2d(x, der)
%% RHS FUNCTION FOR NONLINEAR MONOTONE PROBLEM

if (nargin == 1) || (0 == sum(der))
    t1 = pi .* x(:,1);
    t2 = cos(t1);
    t3 = t2 .^ 2;
    t4 = pi .^ 2;
    t5 = t3 .* t4;
    t6 = (x(:,2) .^ 2);
    t7 = 1 - x(:,2);
    t8 = t7 .^ 2;
    t9 = t6 .* t8;
    t12 = sin(t1);
    t15 = t12 .* t7 - t12 .* x(:,2);
    t17 = 0.1e1 + 0.64e2 .* t5 .* t9 + 0.64e2 .* t15 .^ 2;
    t18 = t17 .^ (-0.125e1);
    t24 = t2 .* pi;
    t37 = t17 .^ (-0.25e0);
    t39 = (0.1e1 + t37) .* t12;
    out = 0.200e1 .* t18 .* (-0.128e3 .* t2 .* t4 .* pi .* t9 .* t12 + 0.128e3 .* t15 .* (t24 .* t7 - t24 .* x(:,2))) .* t2 .* pi .* x(:,2) .* t7 + 0.8e1 .* t39 .* t4 .* x(:,2) .* t7 + 0.200e1 .* t18 .* (0.128e3 .* t5 .* x(:,2) .* t8 - 0.128e3 .* t5 .* t6 .* t7 - 0.256e3 .* t15 .* t12) .* t15 + 0.16e2 .* t39;

else 
    out = zeros(size(x,1), 1);
end

end