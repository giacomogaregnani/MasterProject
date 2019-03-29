function [ out ] = operator_mon3d(x,xi,k)
%% OPERATOR FOR NONLINEAR MONOTONE PROBLEM

t1 = (xi(:,1) .^ 2);
t2 = (xi(:,2) .^ 2);
t3 = (xi(:,3) .^ 2);
t5 = (1 + t1 + t2 + t3) .^ (-0.25e0);

if k == 1
    out = (1 + t5) .* xi(:,1);
    
elseif k == 2
    out = (1 + t5) .* xi(:,2);
    
elseif k == 3
    out = (1 + t5) .* xi(:,3);
    
else
    error('wrong argument');
end

end