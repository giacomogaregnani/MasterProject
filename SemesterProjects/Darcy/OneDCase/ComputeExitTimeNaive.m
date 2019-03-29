function [ExpTau, ExpPhi, t] = ComputeExitTimeNaive(X0,f,g,Bounds,BoundCond,W,Time)
% Computes the EXIT TIME and EXIT PROBABILITY using the DEM method

tic
if BoundCond(2) == 0
    if X0 >= Bounds(2) || X0 <= Bounds(1)
        ExpTau = 0;
        ExpPhi = 1;
        return
    end
elseif BoundCond(2) == 1
    if X0 <= Bounds(1)
        ExpTau = 0;
        ExpPhi = 1;
        return
    end
end

[M,N] = size(W);
h = (Time(2)-Time(1))/(N-1);
tau = Time(2) * ones(M,1);
phi = zeros(M,1);

if BoundCond(2) == 0
    for j = 1:M
        w = W(j,:);
        x = X0;
        for i = 2:N
            x = EMOneStep(x,f,g,w(i)-w(i-1),h);
            if x >= Bounds(2) || x <= Bounds(1)
                tau(j) = h*(i-1);
                phi(j) = 1;
                break
            end
        end
    end
    
elseif BoundCond(2) == 1
    for j = 1:M
        w = W(j,:);
        x = X0;
        for i = 2:N
            x = EMOneStep(x,f,g,w(i)-w(i-1),h);
            if x <= Bounds(1)
                tau(j) = h*(i-1);
                phi(j) = 1;
                break
            elseif x >= Bounds(2)
                x = 2*Bounds(2) - x;
            end
        end
    end
end

ExpTau = mean(tau);
ExpPhi = mean(phi);
t = toc;

end