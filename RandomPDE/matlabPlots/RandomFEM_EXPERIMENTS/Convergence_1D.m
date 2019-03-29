clc; clear; close all
%%

% f1 = @(x) x .^ 0;
f1 = @(x) x.^3 .* (x >= 0.5);
% f1 = @(x) sin(2 * pi * x);
% f1 = @(x) (x-0.5) .* (x >= 0.5);

% f2 = @(x) (x < 0.25) + (0.5 - x) .* (x >= 0.25) .* (x < 0.5) ...
%     + (x - 0.5) .* (x >= 0.5) .* (x < 0.75) + (x >= 0.75);
f2 = @(x) 0 .* x;
% uEx = @(x) 1 / (4*pi^2) * sin(2 * pi * x);
% grad_uEx = @(x) 1 / (2*pi) * cos(2 * pi * x);
field = @(x) x.^0;

pVec = 1:3;

NVec = floor(2.^(2:1:7));
hVec = 1 ./ NVec;

M = 200;

errL2Mean = cell(length(pVec), 1);
errH1Mean = errL2Mean;


for k = 1 : length(pVec)
    
    errL2Mean{k} = zeros(1, length(NVec));
    errH1Mean{k} = zeros(1, length(NVec));    
    errL2loc = zeros(M, length(hVec));
    errH1loc = zeros(M, length(hVec));
    errL2Interploc = zeros(M, length(hVec));
    errH1Interploc = zeros(M, length(hVec));
    errPointloc = zeros(M, length(hVec));
    
    p = pVec(k);
    
    for i = 1 : length(hVec)
        
        N = NVec(i); h = hVec(i);
        x = linspace(0, 1, N + 1);
        
        %         F = assembleRHS(f, x);
        F = assembleRHS_Hminus1(f1, f2, x);
        A = assembleMatrix(field, x);
        u = [0; A \ F; 0];
        gradu = computeGradient_1D(x, u');
        
        figure
        hold on
        
        meanU = zeros(size(u));
        meanX = zeros(size(x));
        
        parfor jj = 1 : M
            
            % Build perturbed grid
            X = x;
            
            for kk = 2 : length(x)-1
                hBar = min(x(kk)-x(kk-1), x(kk+1)-x(kk));
                X(kk) = x(kk) + 0.1 * hBar^p * (rand - 0.5);
            end
            xTot = sort([X, x(2:end-1)]);
            
            meanX = meanX + X;
            
            % FEM solution on perturbed grid
            %             FTilde = assembleRHS(f, X);
            FTilde = assembleRHS_Hminus1(f1, f2, X);
            ATilde = assembleMatrix(field, X);
            U = [0; ATilde\FTilde; 0];
            gradU = computeGradient_1D(X, U');
            
            % Interpolation of FEM solution on perturbed grid
            UI = interp1(x, u, X)';
            grad_UI = computeGradient_1D(X, UI');
            
            % Exact error computations
            errH1loc(jj, i) = 0; errL2loc(jj, i) = 0;
            errL2Interploc(jj, i) = 0; errH1Interploc(jj, i) = 0;
            eL2 = @(xx) (interp1(x,u,xx) - interp1(X,U,xx)).^2;
            eL2I = @(xx) (interp1(x,u,xx) - interp1(X,UI,xx)).^2;
            for kk = 1 : length(xTot) - 1
                xm = mean([xTot(kk), xTot(kk+1)]);
                errH1loc(jj, i) = errH1loc(jj, i) + (xTot(kk+1) - xTot(kk)) * (gradFunction(x, gradu, xm) - gradFunction(X, gradU, xm))^2;
                errL2loc(jj, i) = errL2loc(jj, i) + h/6 * (eL2(xTot(kk)) + 4*eL2(xm) + eL2(xTot(kk+1)));
                errH1Interploc(jj, i) = errH1Interploc(jj, i) + (xTot(kk+1) - xTot(kk)) * (gradFunction(x, gradU, xm) - gradFunction(X, grad_UI, xm))^2;
                errL2Interploc(jj, i) = errL2Interploc(jj, i) + h/6 * (eL2I(xTot(kk)) + 4*eL2I(xm) + eL2I(xTot(kk+1)));
            end
            
            xm = (x(1:end-1) + x(2:end)) / 2;
            Xm = (X(1:end-1) + X(2:end)) / 2;
            errGradPoints = gradFunction(x, gradu, X) - gradFunction(X, gradU, x);
            errPointLoc(jj, i) = max(abs(errGradPoints));
            
            % Plot
%             if ~mod(jj, 10)
%                 plot(X, U, 'color', 0.8 * ones(3, 1))
%             end    
            
            meanU = meanU + U;            
        end
        
        meanU = meanU / M;
        meanX = meanX / M;
        gradMeanU = computeGradient_1D(meanX, meanU');
        plot(x, u, 'k', 'linewidth', 1)
        plot(meanX, meanU, 'r', 'linewidth', 1)
        
        % Compute mean errors
        eL2MeanFunc = @(xx) (interp1(x, u, xx) - interp1(meanX, meanU, xx)).^2;
        xTot = sort([meanX, x(2:end-1)]);
        for kk = 1 : length(xTot) - 1
            xm = mean([xTot(kk), xTot(kk+1)]);
            errL2Mean{k}(i) = errL2Mean{k}(i) + h/6 * (eL2MeanFunc(xTot(kk)) + 4*eL2MeanFunc(xm) + eL2MeanFunc(xTot(kk+1)));
            errH1Mean{k}(i) = errH1Mean{k}(i) + (xTot(kk+1) - xTot(kk)) * (gradFunction(x, gradu, xm) - gradFunction(meanX, gradMeanU, xm))^2;
        end
        
        disp(['iteration ' num2str(k) '.' num2str(i) '/' num2str(length(pVec)) '.' num2str(length(hVec))])
        
    end
    
    if M > 1
        errL2{k} = mean(sqrt(errL2loc));
        errH1{k} = mean(sqrt(errH1loc));
        errL2Interp{k} = mean(sqrt(errL2Interploc));
        errH1Interp{k} = mean(sqrt(errH1Interploc));
        errPoints{k} = mean(errPointLoc);
    else
        errL2{k} = sqrt(errL2loc);
        errH1{k} = sqrt(errH1loc);
        errL2Interp{k} = sqrt(errL2Interploc);
        errH1Interp{k} = sqrt(errH1Interploc);
        errPoints{k} = errPointLoc;
    end
    
end

%%

enhanced = 1;
fontsizeLAB = getLatexTextSize('normalsize', 'enhanced', enhanced);
fontsizeTICK = getLatexTextSize('small', 'enhanced', enhanced);

W = 6.7; H = 6.7;

markers = 'os<';
linestyles = {'--', '-', ':'};

fig = createFigure(W, H, 'enhanced',enhanced);
for k = 1 : length(pVec)
    l(2*k - 1) = loglog(hVec, errL2{k}, 'k-', 'marker', markers(k), 'markersize', 4);
    hold on
    leg{2*k -1 } = ['$p = ' num2str(pVec(k)) '$'];
    loglog(hVec, errL2Mean{k}, 'r-', 'marker', markers(k), 'markersize', 4);
    l(2*k) = loglog(hVec, 2e-2 * hVec.^(pVec(k)+1) , 'k', 'linestyle', linestyles{k});
    leg{2*k} = ['slope ' num2str(pVec(k)+1)];
end
% legH1 = legend(l, leg, 'position', [0.5697,0.1567,0.3298,0.3344], 'interpreter', 'latex');

set(gca, 'xtick', [1e-2, 1e-1])
xlabel('$h$', 'interpreter', 'latex')
ylabel('$\|u_h - \tilde u_h\|_{L^2}$', 'interpreter', 'latex')
axPos = get(gca, 'innerposition');
% ylim([1e-10, 1e-2])
% export_fig(fig, '../../../Reports/DraftPDE_18/VERSION7/Figures/Convergence1D_L2.eps', '-nocrop')

fig = createFigure(W, H, 'enhanced',enhanced);
for k = 1 : length(pVec)
    l(2*k - 1) = loglog(hVec, errH1{k}, 'k-', 'marker', markers(k), 'markersize', 4);
    hold on
    leg{2*k -1} = ['$p = ' num2str(pVec(k)) '$'];
    loglog(hVec, errH1Mean{k}, 'r-', 'marker', markers(k), 'markersize', 4);
    l(2*k) = loglog(hVec, 0.15 * hVec.^(pVec(k)/2 + 1/2), 'k', 'linestyle', linestyles{k});
    leg{2*k} = ['slope ' num2str((pVec(k)+1)/2)];
end
% legL2 = legend(l, leg, 'position', [0.5439,0.1598,0.3506,0.3360], 'interpreter', 'latex');
set(gca, 'xtick', [1e-2, 1e-1])
xlabel('$h$', 'interpreter', 'latex')
ylabel('$|u_h - \tilde u_h|_{H^1}$', 'interpreter', 'latex')
set(gca, 'innerposition', axPos)
% ylim([5e-6, 5e-2])
% export_fig(fig, '../../../Reports/DraftPDE_18/VERSION7/Figures/Convergence1D_H1.eps', '-nocrop')
