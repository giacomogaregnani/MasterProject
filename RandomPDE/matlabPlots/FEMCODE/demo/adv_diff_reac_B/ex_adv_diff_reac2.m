function xx = ex_adv_diff_reac2(x,d)
%EX_ADV_DIFF_REAC2 Summary of this function goes here
%   Detailed explanation goes here

if nargin<2 || (sum(d)==0)
    % u(x)
    xx =sin(pi*x(:,1)) .* sin(pi*x(:,2));
else
    
    % convert 1d input into 2d input
    if length(d) == 1
       if d == 1
           d = [1 0];
       elseif d == 2
           d = [0 1];
       else
           error('wrong argument of d in ex_adv_diff_reac1');   
       end
    elseif length(d) > 2
        error('wrong dimension of d in ex_adv_diff_reac1');
    end
    
    
    if sum(d) == 1 && d(1) == 1
        % D_1 u(x)
        xx = pi * cos(pi*x(:,1)) .* sin(pi*x(:,2));
    elseif sum(d) == 1 && d(2) == 1
        % D_2 u(x)
        xx = pi * sin(pi*x(:,1)) .* cos(pi*x(:,2));
    elseif sum(d) == 2 && d(1) == 2
        % D_11 u(x)
        xx = -pi^2*sin(pi*x(:,1)).*sin(pi*x(:,2));
    elseif sum(d) == 2 && d(2) == 2
        % D_22 u(x)
        xx = -pi^2*sin(pi*x(:,1)).*sin(pi*x(:,2));
    elseif sum(d) == 2
        % D_12 u(x)
        xx = pi^2*cos(pi*x(:,1)).*cos(pi*x(:,2));
    end

end
    
end
