clc; clear; close all
%%

uEx = [];
isFuncEx = true;
field = @(x) x.^0;

f = @(x) sin(2 * pi * x);
uEx = @(x) 1 / (4 * pi^2) * sin(2 * pi * x);
graduEx = @(x) 1 / (2 * pi) * cos(2 * pi * x);
% f = @(x) x.^2 .* (x >= 0.5);
% field = @(x) 1  + (x>=0.2) .* (x<=0.4) + (x>=0.6) .* (x <= 0.8);

% uEx = @(x) -sin(12*pi*x) .* exp(-(x-0.5).^2 / 0.01);
% graduEx = @(x) exp(-100*(x - 1/2).^2).*sin(12*pi*x).*(200*x - 100) - 12*pi*exp(-100*(x - 1/2).^2).*cos(12*pi*x);
% f = @(x) exp(-100*(x - 1/2).^2).*sin(12*pi*x).*(200*x - 100).^2 - 200*exp(-100*(x - 1/2).^2).*sin(12*pi*x) - ...
%     144*pi^2*exp(-100*(x - 1/2).^2).*sin(12*pi*x) - 24*pi*exp(-100*(x - 1/2).^2).*cos(12*pi*x).*(200*x - 100);

if isempty(uEx)
    
    isFuncEx = false;
    
    disp('Computing reference solution...');
    
    xEx = linspace(0, 1, 5e3+1);
    F = assembleRHS(f, xEx);
    A = assembleMatrix(field, xEx);
    
    uExVec = [0; A \ F; 0];
    uEx = @(xx) interp1(xEx, uExVec, xx);
    gradVec = computeGradient_1D(xEx, uExVec);
    graduEx = @(xx) gradFunction(xEx, gradVec, xx);
    
    disp('Computed reference solution');
    
end

%%
p = 1;

NVec = floor(logspace(0.6, 2.5, 8));
hVec = 1 ./ NVec;

M = 5;

errL2Ex = zeros(1, length(NVec));
errH1Ex = zeros(1, length(NVec));
errL2Loc = zeros(M, length(hVec));
errH1Loc = zeros(M, length(hVec));
errL2SumLoc = zeros(M, length(hVec));
errH1SumLoc = zeros(M, length(hVec));

for i = 1 : length(hVec)
    
    N = NVec(i); h = hVec(i);
    x = linspace(0, 1, N+1);
    
    F = assembleRHS(f, x);
    A = assembleMatrix(field, x);
    u = [0; A \ F; 0];
    gradu = computeGradient_1D(x, u');
    
    for jj = 1 : M
        
        % Build perturbed grid
        X = x;
        
        for kk = 2 : length(x)-1
            hBar = min(x(kk)-x(kk-1), x(kk+1)-x(kk));
            X(kk) = x(kk) + hBar^p * (rand - 0.5);
        end
        xTot = sort([X, x(2:end-1)]);
        
        % FEM solution on perturbed grid
        FTilde = assembleRHS(f, X);
        ATilde = assembleMatrix(field, X);
        U = [0; ATilde\FTilde; 0];
        gradU = computeGradient_1D(X, U');
        
        % FEM solution on sum space
        FSum = assembleRHS(f, xTot);
        ASum = assembleMatrix(field, xTot);
        USum = [0; ASum\FSum; 0];
        gradUSum = computeGradient_1D(xTot, USum');
%         errL2SumLoc(jj, i) = integral(@(xx) (uEx(xx) - interp1(xTot,USum,xx)).^2, 0, 1);
%         errH1SumLoc(jj, i) = integral(@(xx) (graduEx(xx) - gradFunction(xTot, gradUSum, xx)).^2, 0, 1);
        
        % Exact error computations
        eL2 = @(xx) (interp1(x,u,xx) - interp1(X,U,xx)).^2;
        eL2Sum = @(xx) (interp1(x,u,xx) - interp1(xTot,USum,xx)).^2;
        for kk = 1 : length(xTot) - 1
            xm = mean([xTot(kk), xTot(kk+1)]);
            errH1Loc(jj, i) = errH1Loc(jj, i) + (xTot(kk+1) - xTot(kk))   * (gradFunction(x, gradu, xm) - gradFunction(X, gradU, xm))^2;
            errL2Loc(jj, i) = errL2Loc(jj, i) + (xTot(kk+1) - xTot(kk))/6 * (eL2(xTot(kk)) + 4*eL2(xm) + eL2(xTot(kk+1)));
            errH1SumLoc(jj, i) = errH1SumLoc(jj, i) + (xTot(kk+1) - xTot(kk))   * (gradFunction(x, gradu, xm) - gradFunction(xTot, gradUSum, xm))^2;
            errL2SumLoc(jj, i) = errL2SumLoc(jj, i) + (xTot(kk+1) - xTot(kk))/6 * (eL2Sum(xTot(kk)) + 4*eL2Sum(xm) + eL2Sum(xTot(kk+1)));
        end
        
    end
    
    % Compute exact errors
    if isFuncEx
        errL2Ex(i)  = sqrt(integral(@(xx) (uEx(xx) - interp1(x,u,xx)).^2, 0, 1));
        errH1Ex(i)  = sqrt(integral(@(xx) (graduEx(xx) - gradFunction(x, gradu, xx)).^2, 0, 1));
    else
        eL2 = @(xx) (interp1(x,u,xx) - uEx(xx)).^2;
        for kk = 1 : length(xEx) - 1
            xm = mean([xEx(kk), xEx(kk+1)]);
            errH1Ex(i) = errH1Ex(i) + (xEx(kk+1) - xEx(kk))   * (gradFunction(x, gradu, xm) - graduEx(xm))^2;
            errL2Ex(i) = errL2Ex(i) + (xEx(kk+1) - xEx(kk))/6 * (eL2(xEx(kk)) + 4*eL2(xm) + eL2(xEx(kk+1))) ;
        end
        errH1Ex(i) = sqrt(errH1Ex(i));
        errL2Ex(i) = sqrt(errL2Ex(i));
    end
    
    disp(['iteration ' num2str(i) '/' num2str(length(hVec))])
    
    if M > 1
        errL2 = sqrt(mean(errL2Loc));
        errH1 = sqrt(mean(errH1Loc));
        errL2Sum = sqrt(mean(errL2SumLoc));
        errH1Sum = sqrt(mean(errH1SumLoc));
    else
        errL2 = sqrt(errL2Loc);
        errH1 = sqrt(errH1Loc);
        errL2Sum = sqrt(errL2SumLoc);
        errH1Sum = sqrt(errH1SumLoc);
    end
end


%%

enhanced = 2;
fontsizeLAB = getLatexTextSize('normalsize', 'enhanced', enhanced);
fontsizeTICK = getLatexTextSize('footnotesize', 'enhanced', enhanced);

W = 6.7; H = 6.7;

% L2 error
fig = createFigure(W, H, 'enhanced',enhanced);
loglog(hVec, errL2, 'k<-', 'markersize', 4);
hold on
loglog(hVec, errL2Ex, 'ks-', 'markersize', 4);
loglog(hVec, errL2Sum, 'kd-', 'markersize', 4);
loglog(hVec, 10 * hVec.^2, 'k--')
legend({'$E \|u_h - \tilde u_h\|_{L^2}$', '$\|u_h - u\|_{L^2}$', '$E \|u_h - u_h^+\|_{L^2}$', 'slope 2'}, 'location', 'se', 'interpreter', 'latex', 'fontsize', fontsizeTICK)
% legend({'$E \|u_h - \tilde u_h\|_{L^2}$', '$\|u_h - u\|_{L^2}$', 'slope 2'}, 'location', 'nw', 'interpreter', 'latex')
xlabel('$h$', 'interpreter', 'latex')
ylabel('$L^2$ error', 'interpreter', 'latex')
inPos = get(gca, 'innerposition');
%export_fig(fig, '../../../Reports/DraftPDE_18/VERSION7/Figures/ErrEst1D_ConvL2.eps', '-nocrop')

% L2 effectivity index
fig = createFigure(W, H, 'enhanced',enhanced);
semilogx(hVec, errL2 ./ errL2Ex, 'k<-', 'markersize', 4);
xlabel('$h$', 'interpreter', 'latex')
ylabel('$L^2$ effectivity index', 'interpreter', 'latex')
set(gca, 'innerposition', inPos)
%export_fig(fig, '../../../Reports/DraftPDE_18/VERSION7/Figures/ErrEst1D_EffectL2.eps', '-nocrop')

% H1 error
fig = createFigure(W, H, 'enhanced',enhanced);
loglog(hVec, errH1, 'k<-', 'markersize', 4);
hold on
loglog(hVec, errH1Ex, 'ks-', 'markersize', 4);
loglog(hVec, errH1Sum, 'kd-', 'markersize', 4);
loglog(hVec, 50 * hVec, 'k--')
legend({'$E|u_h - \tilde u_h|_{H^1}$', '$|u_h - u|_{H^1}$', '$E |u_h - u_h^+|_{H^1}$', 'slope 1'}, 'location', 'se', 'interpreter', 'latex', 'fontsize', fontsizeTICK)
% legend({'$E|u_h - \tilde u_h|_{H^1}$', '$|u_h - u|_{H^1}$', 'slope 1'}, 'location', 'nw', 'interpreter', 'latex')
set(gca, 'innerposition', inPos)
xlabel('$h$', 'interpreter', 'latex')
ylabel('$H^1$ error', 'interpreter', 'latex')
%export_fig(fig, '../../../Reports/DraftPDE_18/VERSION7/Figures/ErrEst1D_ConvH1.eps', '-nocrop')

% H1 effectivity index
fig = createFigure(W, H, 'enhanced',enhanced);
set(gca, 'innerposition', inPos)
semilogx(hVec, errH1 ./ errH1Ex, 'k<-', 'markersize', 4);
xlabel('$h$', 'interpreter', 'latex')
ylabel('$H^1$ effectivity index', 'interpreter', 'latex')
%export_fig(fig, '../../../Reports/DraftPDE_18/VERSION7/Figures/ErrEst1D_EffectH1.eps', '-nocrop')
