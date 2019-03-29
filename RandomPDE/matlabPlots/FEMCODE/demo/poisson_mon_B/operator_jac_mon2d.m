function [ out ] = operator_jac_mon2d(x,xi,k,l)
%% JACOBIAN MATRIX OF OPERATOR FOR NONLINEAR MONOTONE PROBLEM

t1 = (xi(:,1) .^ 2);
t2 = (xi(:,2) .^ 2);
t3 = 1 + t1 + t2;
t4 = t3 .^ (-0.125e1);
t7 = t3 .^ (-0.25e0);


if k == 1 && l == 1
    out = -0.50e0 .* t4 .* t1 + 0.1e1 + t7;

elseif k == 1 && l == 2
    out = -0.50e0 .* t4 .* xi(:,2) .* xi(:,1);
    
elseif k == 2 && l == 1
    out = -0.50e0 .* t4 .* xi(:,2) .* xi(:,1);
    
elseif k == 2 && l == 2
    out = -0.50e0 .* t4 .* t2 + 0.1e1 + t7;
    
else
    out = zeros(size(x,1),1);
end

end