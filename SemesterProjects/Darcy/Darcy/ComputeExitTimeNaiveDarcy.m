function [ExpTau,ExpPhi,t] = ComputeExitTimeNaiveDarcy(X0,g,Bounds,BoundCond,W,Time,Ux,Uy,delta)
% Performs DEM in DARCY setting

tic
sigma = det(g(0,0));

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
        x = X0;
        for i = 2:N
            % find where I am and find the value of the velocity field
            index = [ceil((x(1)+1)/delta),ceil((x(2)+1)/delta)];
            u = [Ux(index(1),index(2)); Uy(index(1),index(2))];
            x = EMOneStepDarcy(x,u,sigma,w(:,i)-w(:,i-1),h);
            if x(1) >= Bounds(1,2) || x(1) <= Bounds(1,1) || x(2) >= Bounds(2,2) || x(2) <= Bounds(2,1)
                tau(j) = h*(i-1);
                phi(j) = 1;
                break
            end
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
        x = X0;
        for i = 2:N
            % find where I am and find the value of the velocity field
            index = [ceil((x(1)+1)/delta),ceil((x(2)+1)/delta)];
            u = [Ux(index(1),index(2)); Uy(index(1),index(2))];
            x = EMOneStepDarcy(x,u,sigma,w(:,i)-w(:,i-1),h);
            if x(1) >= Bounds(1,2) || x(1) <= Bounds(1,1)
                tau(j) = h*(i-1);
                phi(j) = 1;
                break
            elseif x(2) < Bounds(2,1)
                x(2) = 2*Bounds(2,1) - x(2);
            elseif x(2) > Bounds(2,2)
                x(2) = 2*Bounds(2,2) - x(2);
            end
        end
    end
end
ExpTau = mean(tau);
ExpPhi = mean(phi);

t = toc;
end