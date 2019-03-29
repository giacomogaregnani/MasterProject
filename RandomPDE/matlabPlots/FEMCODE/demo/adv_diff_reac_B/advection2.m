function xx = advection2(x,d)
%ADVECTION2 Summary of this function goes here
%   Detailed explanation goes here

if d==1
    xx = -x(:,1);
elseif d==2
    xx = 1./(x(:,1).^2 + 0.1);
else
    error('wrong arguments for velocity_field1')
end

end

