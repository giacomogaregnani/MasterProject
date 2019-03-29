function xx = ex_adv_diff_reac1_3D(x,d)
%% EXACT SOLUTION FOR A 3D ELLIPTIC PROBLEM

if nargin<2 || (sum(d)==0)
    % u(x)
    xx =sin(pi*x(:,1)) .* sin(pi*x(:,2)) .*  sin(pi*x(:,3));
else
    
    % convert 1d input into 2d input
    if length(d) == 1
       if d == 1
           d = [1 0 0];
       elseif d == 2
           d = [0 1 0];
       elseif d == 3
           d = [0 0 1];
       else
           error('wrong argument of d in ex_adv_diff_reac1_3D');   
       end
    elseif length(d) > 3
        error('wrong dimension of d in ex_adv_diff_reac1_3D');
    end
    
    
    if sum(d) == 1 && d(1) == 1
        % D_1 u(x)
        xx = pi * cos(pi*x(:,1)) .* sin(pi*x(:,2)) .* sin(pi*x(:,3));
    elseif sum(d) == 1 && d(2) == 1
        % D_2 u(x)
        xx = pi * sin(pi*x(:,1)) .* cos(pi*x(:,2)) .* sin(pi*x(:,3));
    elseif sum(d) == 1 && d(3) == 1
        % D_3 u(x)
        xx = pi * sin(pi*x(:,1)) .* sin(pi*x(:,2)) .* cos(pi*x(:,3));
    elseif sum(d) == 2 && d(1) == 2
        % D_11 u(x)
        xx = -pi^2*sin(pi*x(:,1)).*sin(pi*x(:,2)).* sin(pi*x(:,3));
    elseif sum(d) == 2 && d(2) == 2
        % D_22 u(x)
        xx = -pi^2*sin(pi*x(:,1)).*sin(pi*x(:,2)) .* sin(pi*x(:,3));
    elseif sum(d) == 2 && d(3) == 2
        % D_33 u(x)
        xx = -pi^2*sin(pi*x(:,1)).*sin(pi*x(:,2)) .* sin(pi*x(:,3));
    else
        error('this derivative is not implemented yet');
    end

end
    
end

