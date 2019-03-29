function [ out ] = tensor_nonmon3d(x,s,k,l)
%% TENSOR FOR NONLINEAR NONMONOTONE PROBLEM

if k == 1 && l == 1
    out = 1+x(:,1).*sin(pi*s);
    
elseif k == 2 && l == 2
    out = 2+atan(s);

elseif k==3 && l ==3
    out = 1+x(:,3).*cos(pi*s);

else
    out = zeros(size(x,1),1);
end

end
