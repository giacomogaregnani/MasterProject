clc; clear; close all
%%

uEx = [];
kappa = @(x) x.^0;
% kappa = @(x) 1 + 2.3 * (x >= 0.3) .* (x <= 0.7);

% uEx = @(x) 1 / (4 * pi^2) * sin(2 * pi * x);
% graduEx = @(x) 1 / (2 * pi) * cos(2 * pi * x);
% f = @(x) sin(2 * pi * x);

% uEx = @(x) -sin(4*pi*x) .* (cos(2*pi*x) + sin(4*pi*x));
% graduEx = @(x) -sin(4*pi*x).*(4*pi*cos(4*pi*x) - 2*pi*sin(2*pi*x)) - 4*pi*cos(4*pi*x).*(cos(2*pi*x) + sin(4*pi*x));
% f = @(x) 8*pi*cos(4*pi*x).*(4*pi*cos(4*pi*x) - 2*pi*sin(2*pi*x)) - sin(4*pi*x).*(4*pi^2*cos(2*pi*x) ...
%     + 16*pi^2*sin(4*pi*x)) - 16*pi^2*sin(4*pi*x).*(cos(2*pi*x) + sin(4*pi*x));

% uEx = @(x) -sin(2*pi*x) .* exp(-100 * (x-0.5).^2);
% graduEx = @(x) exp(-100*(x - 1/2).^2).*sin(2*pi*x).*(200*x - 100) - 2*pi*exp(-100*(x - 1/2).^2).*cos(2*pi*x);
% f = @(x) exp(-100*(x - 1/2).^2).*sin(2*pi*x).*(200*x - 100).^2 - 200*exp(-100*(x - 1/2).^2).*sin(2*pi*x) - ...
%     4*pi^2*exp(-100*(x - 1/2).^2).*sin(2*pi*x) - 4*pi*exp(-100*(x - 1/2).^2).*cos(2*pi*x).*(200*x - 100);

uEx = @(x) -sin(24*pi*x) .* exp(-(x-0.5).^2 / 0.01);
graduEx = @(x) exp(-100*(x - 1/2).^2).*sin(24*pi*x).*(200*x - 100) - 24*pi*exp(-100*(x - 1/2).^2).*cos(24*pi*x);
f = @(x) exp(-100*(x - 1/2).^2).*sin(24*pi*x).*(200*x - 100).^2 - 200*exp(-100*(x - 1/2).^2).*sin(24*pi*x) - ...
    576*pi^2*exp(-100*(x - 1/2).^2).*sin(24*pi*x) - 48*pi*exp(-100*(x - 1/2).^2).*cos(24*pi*x).*(200*x - 100);

% f = @(x) x.^3 .* (x >= 0.5);
% f2 = @(x) x.^3 .* (x >= 0.5);
% f1 = @(x) (x < 0.25) + (0.5 - x) .* (x >= 0.25) .* (x < 0.5) ...
%     + (x - 0.5) .* (x >= 0.5) .* (x < 0.75) + (x >= 0.75);
% f2 = @(x) 0 * x;

xFine = 0:1e-4:1;

if isempty(uEx)
    AEx = assembleMatrix(kappa, xFine);
    FEx = assembleRHS(f, xFine);
    uExVec = [0; AEx \ FEx; 0];
    uGrad = computeGradient_1D(xFine, uExVec);
    uEx = @(xx) interp1(xFine, uExVec, xx);
    graduEx = @(xx) gradFunction(xFine, uGrad, xx);
end

solNorm = sqrt(integral(@(x) uEx(x).^2, 0, 1));

%%

quad = false;
errTrue = []; errEst = []; N = [];
x = linspace(0, 1, 20);
flag = true;
p = 1;
M = 5e2;
it = 0;
tol = 1e-3;

while flag
    
    it = it + 1;
    
    N(it) = length(x)-1;
    h = max(x(2:end) - x(1:end-1));
    
    err = zeros(1, N(it));
    
    % Deterministic solution
    if quad
        F = assembleRHS_quadrature(f, x);
    else
        F = assembleRHS(f, x);
    end
    A = assembleMatrix(kappa, x);
    u = [0; A \ F; 0];
    for kk = 1 : N(it)
        err(kk) = sqrt(integral(@(xx) (uEx(xx) - interp1(x, u, xx)).^2, x(kk), x(kk+1)));
    end
    errTrue(it) = sqrt(sum(err.^2)) / solNorm;
    
    % Probabilistic estimation
    errProbSqd = zeros(M, N(it));
    %     figure
    %     hold on
    parfor jj = 1 : M
        X = x;
        for kk = 2 : length(x)-1
            hBar = min(x(kk)-x(kk-1), x(kk+1)-x(kk));
            X(kk) = x(kk) + hBar^p * (rand - 0.5);
        end
        if quad
            FTilde = assembleRHS_quadrature(f, X);
        else
            FTilde = assembleRHS(f, X);
        end
        ATilde = assembleMatrix(kappa, X);
        uTilde = [0; ATilde \ FTilde; 0];
        tmp = zeros(1, N(it));
        for kk = 1 : N(it)
            tmp(kk) = integral(@(xx) (interp1(x, u, xx) - interp1(X, uTilde, xx)).^2, X(kk), X(kk+1));
        end
        errProbSqd(jj, :) = tmp;
        %         plot(X, uTilde, 'color', 0.7 * ones(1, 3));
    end
    %     plot(x, u, 'k--', 'linewidth', 2)
    %     plot(xFine, uEx(xFine), 'k', 'linewidth', 2)
    if M > 1
        errProbSqd = mean(errProbSqd);
    end
    errEst(it) = sqrt(sum(errProbSqd)) / solNorm;
    errProb = sqrt(errProbSqd);
    
    % Refine the grid       
    if errEst(it) < tol
        flag = false;
    else
        [xNew, nNewPoints, nOldPoints]  = refineGrid(x, errProb, tol, solNorm);
    end
        
    % Plot
    figure
    semilogy((x(2:end) + x(1:end-1))/2, err, 'k')
    hold on
    semilogy((x(2:end) + x(1:end-1))/2, errProb, 'b')
    legend({'true', 'est'}, 'interpreter', 'latex')
    
    figure
    plot(x, zeros(1, N(it)+1), 'r.')
    
    x = xNew;
    
    display(['iteration ', num2str(it), ', num elem ', num2str(N(it)), ', error ', num2str(errEst(it)) ...
        ', new points ', num2str(nNewPoints), ', old points ' , num2str(nOldPoints)]);
    
end

loglog(N, errEst, 'ks-')
hold on
loglog(N, errTrue, 'k<-')
loglog(N, tol * ones(size(N)), 'k--');
legend({'estimate', 'true', 'tol'}, 'location', 'ne', 'interpreter', 'latex')



