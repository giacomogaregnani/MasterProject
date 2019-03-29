function [ out ] = tensor_jac_nonmon3d(x,s,k,l)
%% JACOBIAN MATRIX OF TENSOR FOR NONLINEAR NONMONOTONE PROBLEM

if k == 1 && l == 1
    out = pi*x(:,1).*cos(pi*s);
    
elseif k == 2 && l == 2
    out = 1./(1+s.^2);

elseif k == 3 && l == 3
    out = -pi*x(:,3).*sin(pi*s);

else
    out = zeros(size(x,1),1);
end

end
