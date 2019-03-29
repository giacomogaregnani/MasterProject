function xx = diffusion2(x,k,l,d)
%DIFFUSION2 Summary of this function goes here
%   Detailed explanation goes here

if nargin == 3 || sum(d) == 0
    
    if k==1 && l==1
        % a11
        xx = ones(size(x,1),1);
        return;
    elseif k==1 && l==2
        % a12
        xx = zeros(size(x,1),1);
        return;
    elseif k==2 && l==1
        % a21
        xx = zeros(size(x,1),1);
        return;
    elseif k==2 && l==2
        % a22
        xx = ones(size(x,1),1);
        return
    end
    
elseif sum(d) == 1 && d(1) == 1
    
    if k==1 && l==1
        % D_1 a11(x)
        xx = zeros(size(x,1),1);
        return;
    elseif k==1 && l==2
        % D_1 a12(x)
        xx = zeros(size(x,1),1);
        return;
    end
    
elseif sum(d) == 1 && d(2) == 1
    
    if k==2 && l==1
        % D_2 a21(x)
        xx = zeros(size(x,1),1);
        return;
    elseif k==2 && l==2
        % D_2 a22(x)
        xx = zeros(size(x,1),1);
        return;
    end
    
else
    error('wrong arguments in diffusion1');
end

end

