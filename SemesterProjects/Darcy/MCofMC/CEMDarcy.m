function [ExpTau, ExpPhi, t] = CEMDarcy(X0, Sigma, Bounds, BoundCond, W, Time, Ux, Uy, delta)
% Performs CEM in DARCY case

tic

if BoundCond == 0
    if X0(1) >= Bounds(1,2) || X0(1) <= Bounds(1,1) || X0(2) >= Bounds(2,2) || X0(2) <= Bounds(2,1)
        ExpTau = 0;
        return
    end
    
    [TwoM,N] = size(W);
    M = TwoM/2;
    h = (Time(2)-Time(1))/(N-1);
    tau = Time(2) * ones(M,1);
    phi = zeros(M,1);
    
    for j = 1:M
        w = W(2*j-1:2*j,:);
        xOld = X0;
        for i = 2:N

            index = [ceil((xOld(1)+1)/delta),ceil((xOld(2)+1)/delta)];
            u = [Ux(index(1),index(2)); Uy(index(1),index(2))];
            xNew = EMOneStepDarcy(xOld,u,Sigma,w(:,i)-w(:,i-1),h);
            
            if xNew(1) >= Bounds(1,2) || xNew(1) <= Bounds(1,1) || xNew(2) >= Bounds(2,2) || xNew(2) <= Bounds(2,1)
                tau(j) = h*(i-1);
                phi(j) = 1;
                break
            else
                p = ComputeExitProbability(xOld,xNew,Sigma,h);
                u = rand(1,1);
                if  isempty(find(p > u,1)) == 0
                    tau(j) = h*(i-1);
                    phi(j) = 1;
                    break
                end
            end
            xOld = xNew;
        end
    end
    
else
    if X0(1) >= Bounds(1,2) || X0(1) <= Bounds(1,1)
        ExpTau = 0;
        return
    end
    
    [TwoM,N] = size(W);
    M = TwoM/2;
    h = (Time(2)-Time(1))/(N-1);
    tau = Time(2) * ones(M,1);
    phi = zeros(M,1);
    
    for j = 1:M
        w = W(2*j-1:2*j,:);
        xOld = X0;
        for i = 2:N

            index = [ceil((xOld(1)+1)/delta),ceil((xOld(2)+1)/delta)];
            u = [Ux(index(1),index(2)); Uy(index(1),index(2))];
            xNew = EMOneStepDarcy(xOld,u,Sigma,w(:,i)-w(:,i-1),h);
            
            if xNew(1) >= Bounds(1,2) || xNew(1) <= Bounds(1,1)
                tau(j) = h*(i-1);
                phi(j) = 1;
                break
            elseif xNew(2) < Bounds(2,1)
                xNew(2) = 2*Bounds(2,1) - xNew(2);
            elseif xNew(2) > Bounds(2,2)
                xNew(2) = 2*Bounds(2,2) - xNew(2);
            else
                p = ComputeExitProbability(xOld,xNew,Sigma,h);
                p = [p(1),0,p(3),0];
                u = rand(1,1);
                if  isempty(find(p > u,1)) == 0
                    tau(j) = h*(i-1);
                    phi(j) = 1;
                    break
                end
            end
            xOld = xNew;
        end
    end
end

t = toc;
ExpTau = mean(tau);
ExpPhi = mean(phi);
end
