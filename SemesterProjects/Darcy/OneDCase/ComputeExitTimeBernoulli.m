function [ExpTau, ExpPhi, t] = ComputeExitTimeBernoulli(X0, f, g, Bounds, BoundCond, W, Time)
% Computes the EXIT TIME and EXIT PROBABILITY using the CEM method

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
sigma = g(1);

if BoundCond(2) == 0
    for j = 1:M
        w = W(j,:);
        xOld = X0;
        for i = 2:N
            xNew = EMOneStep(xOld,f,g,w(i)-w(i-1),h);
            if xNew >= Bounds(2) || xNew <= Bounds(1)
                tau(j) = h*(i-1);
                phi(j) = 1;
                break
            else
                p = TwoSidedExitProbability(xOld,xNew,sigma,h,Bounds);
                unif = rand(1,1);
                if unif <= p
                    tau(j) = h*(i-1);
                    phi(j) = 1;
                    break
                end
            end
            xOld = xNew;
        end
    end
    
elseif BoundCond(2) == 1
    for j = 1:M
        w = W(j,:);
        xOld = X0;
        for i = 2:N
            xNew = EMOneStep(xOld,f,g,w(i)-w(i-1),h);
            if xNew <= Bounds(1)
                tau(j) = h*(i-1);
                phi(j) = 1;
                break
            elseif xNew > Bounds(2)
                xNew = 2*Bounds(2) - xNew;
            else
                p = exp(-2 * ((xOld - Bounds(1)) * (xNew - Bounds(1))) / (sigma^2 * h));
                unif = rand(1,1);
                if unif <= p
                    tau(j) = h*(i-1);
                    phi(j) = 1;
                    break
                end
            end
            xOld = xNew;
        end
    end
end

ExpTau = mean(tau);
ExpPhi = mean(phi);
t = toc;

end