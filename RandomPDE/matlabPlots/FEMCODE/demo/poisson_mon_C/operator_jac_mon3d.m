function [ out ] = operator_jac_mon3d(x,xi,k,l)
%% JACOBIAN MATRIX OF OPERATOR FOR NONLINEAR MONOTONE PROBLEM

t1 = (xi(:,1) .^ 2);
t2 = (xi(:,2) .^ 2);
t3 = (xi(:,3) .^ 2);
t4 = 1 + t1 + t2 + t3;
t5 = t4 .^ (-0.125e1);
t8 = t4 .^ (-0.25e0);

if k == 1 && l == 1
    out = -0.50e0 .* t5 .* t1 + 0.1e1 + t8;

elseif k == 1 && l == 2
    out = -0.50e0 .* t5 .* xi(:,1) .* xi(:,2);
    
elseif k == 1 && l == 3
    out = -0.50e0 .* t5 .* xi(:,3) .* xi(:,1);
    
elseif k == 2 && l == 1
    out = -0.50e0 .* t5 .* xi(:,1) .* xi(:,2);
    
elseif k == 2 && l == 2
    out = -0.50e0 .* t5 .* t2 + 0.1e1 + t8;
    
elseif k == 2 && l == 3
    out = -0.50e0 .* t5 .* xi(:,2) .* xi(:,3);
    
elseif k == 3 && l == 1
    out = -0.50e0 .* t5 .* xi(:,3) .* xi(:,1);
    
elseif k == 3 && l == 2
    out = -0.50e0 .* t5 .* xi(:,2) .* xi(:,3);
    
elseif k == 3 && l == 3
    out = -0.50e0 .* t5 .* t3 + 0.1e1 + t8;
    
else
    out = zeros(size(x,1),1);
end

end