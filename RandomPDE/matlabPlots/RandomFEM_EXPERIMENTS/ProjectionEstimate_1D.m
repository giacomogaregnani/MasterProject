clc; clear; close all
%%

% POINTS = zeros(1, 11);
% POINTS(2:2:end-1) = 1;
% POINTS = [0, POINTS, 0];
% X = [-1, linspace(0, 1, length(POINTS)-2), 2];
% f = @(xx) xx + 0.1 * interp1(X, POINTS, xx);

% f = @(x) sin(2*pi*x);
f = @(x) (x - 0.5).^2 .* (x >= 0.5);
% f = @(x) x.^2 .* (x <= 0.47) + (x - 0.94).^2 .* (x > 0.47);
% f = @(x) 1 + 0.01 * x.^(3/2);
% f = @(x) sin(100 * pi * x);
% f = @(x) x .* (x < 0.23) + (0.5 - x) .* (x >= 0.23) .* (x < 0.5) ...
%        + (x - 0.5) .* (x >= 0.5) .* (x < 0.75) + (1 - x) .* (x >= 0.75);
% f = @(x) x.^2 .* (x < 0.5) + (x-1).^2 .* (x >= 0.5);
% f = @(x) exp(-10*(x - 1/2).^2).*sin(24*pi*x).*(20*x - 10).^2 - 20*exp(-10*(x - 1/2).^2).*sin(24*pi*x) - ...
%     576*pi^2*exp(-10*(x - 1/2).^2).*sin(24*pi*x) - 48*pi*exp(-10*(x - 1/2).^2).*cos(24*pi*x).*(20*x - 10);
% f = @(x) interp1([-1 0.1, 0.2, 0.3 2], [0, 0, 1, 0, 0], x); 

NVec = sort(unique([4:10, floor(2.^(2:0.4:7))]));
hVec = 1 ./ NVec;

M = 100;
p = 1;

errLoopH1 = zeros(1, M);
errLoopL2 = zeros(1, M);
errLoopLinf = zeros(1, M);
errGradPoint = zeros(1, M);

for i = 1 : length(hVec)
    
    h = hVec(i);
    x = 0:h:1;
    disp([num2str(i) '/' num2str(length(hVec)) ': h = ' num2str(h) ' N = ' num2str(length(x)) ]);
    
    F = f(x); 
  
    grad_F = computeGradient_1D(x, F);
    grad_Norm = sum(h * grad_F.^2)^(1/2);
    
%     figure
%     hold on
    parfor j = 1 : M 
        
%         vecc = 0.2 * ones(size(x)); vecc(1:2:end) = -vecc(1:2:end);
        X = x + h^p * (rand(size(x)) - 0.5); 
        X(1) = 0; 
        X(end) = 1;
            
        FF = interp1(x,F,X);
        grad_FF = computeGradient_1D(X, FF);
        
        [xTot, indices] = sort([X, x(2:end-1)]);          
        FFSumSpace = [F, FF(2:end-1)]; FFSumSpace = FFSumSpace(indices);
        grad_FFSumSpace = computeGradient_1D(xTot, FFSumSpace);
        
        errLoopH1(j) = 0; errLoopL2(j) = 0;
        eL2 = @(xx) (interp1(x,F,xx) - interp1(X,FF,xx)).^2;
        for k = 1 : length(xTot) - 1
            xm = mean([xTot(k), xTot(k+1)]);
            errLoopH1(j) = errLoopH1(j) + (xTot(k+1) - xTot(k)) * (gradFunction(x, grad_F, xm) - gradFunction(X, grad_FF, xm))^2;
%             errLoopH1(j) = errLoopH1(j) + (xTot(k+1) - xTot(k)) * (gradFunction(x, grad_F, xm) - gradFunction(xTot, grad_FFSumSpace, xm))^2;
            errLoopL2(j) = errLoopL2(j) + h/6 * (eL2(xTot(k)) + 4*eL2(xm) + eL2(xTot(k+1)));
        end
        errLoopLinf(j) = max(abs(interp1(x,F,x) - interp1(X,FF,x)));
                
        xm = (x(1:end-1) + x(2:end)) / 2;
        Xm = (X(1:end-1) + X(2:end)) / 2;
        errGradPoints = gradFunction(x, grad_F, xm) - gradFunction(X, grad_FF, Xm);
        errGradPoint(j) = max(abs(errGradPoints));
       
    end
%     plot(x, F, 'k');
   
    errH1(i) = mean(sqrt(errLoopH1) / grad_Norm);
    errL2(i) = mean(sqrt(errLoopL2));
    errPoint(i) = mean(errGradPoint);
    errLinf(i) = mean(errLoopLinf);
    
end

%%
figure
loglog(hVec, errH1, 'o-');
hold on
loglog(hVec, errPoint, 'o-');
loglog(hVec, errL2, 'o-');
loglog(hVec, errLinf, 'o-');
loglog(hVec, grad_Norm*hVec.^((p-1)/2), 'k-.')
loglog(hVec, grad_Norm*hVec.^((p+1)/2), 'k-');
loglog(hVec, grad_Norm*hVec.^(p+1/2), 'k:');
loglog(hVec, grad_Norm*hVec.^(p+1), 'k--');
loglog(hVec, grad_Norm*hVec.^(p), 'm')
loglog(hVec, grad_Norm*hVec.^(p/2), 'm--')
loglog(hVec, grad_Norm*hVec.^(p-1), 'm:')
legend('H1', 'point', 'L2', 'inf', '(p-1)/2', '(p+1)/2', 'p+1/2', 'p+1', 'p', 'p/2', 'p-1', 'location', 'best')

figure
loglog(hVec, errPoint, 'o-');
hold on
loglog(hVec, 0.5 * grad_Norm * hVec.^p .* abs(log(hVec)))