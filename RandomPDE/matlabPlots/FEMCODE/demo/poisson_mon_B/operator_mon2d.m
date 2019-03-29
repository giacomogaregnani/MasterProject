function [ out ] = operator_mon2d(x,xi,k)
%% OPERATOR FOR NONLINEAR MONOTONE PROBLEM

t1 = (xi(:,1) .^ 2);
t2 = (xi(:,2) .^ 2);
t4 = (1 + t1 + t2) .^ (-0.25e0);

if k == 1
    out = (1 + t4) .* xi(:,1);
    
elseif k == 2
    out = (1 + t4) .* xi(:,2);
else
    error('wrong argument');
end

end