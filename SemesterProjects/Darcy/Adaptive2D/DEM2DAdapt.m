function [ExpTau, ExpPhi, nStep, H, D] = DEM2DAdapt(X0,f,g,Bounds,BoundCond,M,Time,hmin,h0)
% Estimate EXIT TIME and EXIT PROBABILITY using ADAPTIVE METHOD

nStep = 0;

if BoundCond == 0
    if X0(1) >= Bounds(1,2) || X0(1) <= Bounds(1,1) || X0(2) >= Bounds(2,2) || X0(2) <= Bounds(2,1)
        ExpTau = 0;
        return
    end
    
    tau = Time(2) * ones(M,1);
    phi = zeros(M,1);
    sigma = g(1,1);
    sigma = sigma(1,1);
    l = -log2(hmin / h0);
    nStep = phi;
    
    for j = 1:M
        x = X0;
        time = 0;
        H = [];
        D = [];
        while time < Time(2)
            
            % compute distance
            r1 = x(1) + 1;
            r2 = x(2) + 1;
            r3 = 1 - x(1);
            r4 = 1 -x(2);
            d = min([r1, r2, r3, r4]);
            
            D = [D d];
            
            funcD = (d / ((l + 3) * sigma))^2;
            h = max(hmin, min(h0, funcD));
            
            if time + h > Time(2)
                h = Time(2) - time;
            end
            
            H = [H, h];
            nStep(j) = nStep(j) + 1;
            x = EMOneStep(x,f,sigma,h);
            time = time + h;
            
            if x(1) >= Bounds(1,2) || x(1) <= Bounds(1,1) || x(2) >= Bounds(2,2) || x(2) <= Bounds(2,1)
                tau(j) = time;
                phi(j) = 1;
                time = Time(2);
            end
        end
    end
    
elseif BoundCond == 1
    
    if X0(1) >= Bounds(1,2) || X0(1) <= Bounds(1,1)
        ExpTau = 0;
        return
    end
    
    tau = Time(2) * ones(M,1);
    phi = zeros(M,1);
    sigma = g(1,1);
    sigma = sigma(1,1);
    l = -log2(hmin / h0);
    nStep = phi;
    
    for j = 1:M
        x = X0;
        time = 0;
        H = [];
        D = [];
        
        while time < Time(2)
            
            d = ComputeDistance(x);
            funcD = (d / ((l + 3) * sigma))^2;
            h = max(hmin, min(h0, funcD));
            
            if time + h > Time(2)
                h = Time(2) - time;
            end
            
            H = [H, h];
            D = [D, d];
            x = EMOneStep(x,f,sigma,h);
            nStep(j) = nStep(j) + 1;
            time = time + h;
            
            if x(1) >= Bounds(1,2) || x(1) <= Bounds(1,1)
                tau(j) = time;
                phi(j) = 1;
                time = Time(2);
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
nStep = mean(nStep);

end