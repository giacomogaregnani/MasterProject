function xx = ex_WeakAdvDiffMS2D(x,d)
%% EXACT SOLUTION FOR A 2D ELLIPTIC PROBLEM

if nargin<2 || (sum(d)==0)
    % u(x)
    xx =sin(2*pi*x(:,1)) .* sin(2*pi*x(:,2));
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
        xx = 2*pi * cos(2*pi*x(:,1)) .* sin(2*pi*x(:,2));
    elseif sum(d) == 1 && d(2) == 1
        % D_2 u(x)
        xx = 2*pi * sin(2*pi*x(:,1)) .* cos(2*pi*x(:,2));
    else
        error('not implemented');
    end

end
    
end

