clc; clear; close all
%%

% u_Ex = @(x) 1 / (4 * pi^2) * sin(2 * pi * x);
% gradu_Ex = @(x) 1 / (2 * pi) * cos(2 * pi * x);
% f = @(x) sin(2 * pi * x);

enhanced = 1;
fontsizeLAB = getLatexTextSize('normalsize', 'enhanced', enhanced);
fontsizeTICK = getLatexTextSize('small', 'enhanced', enhanced);
W = 6; H = 6;

%%

uEx = @(x) -sin(12*pi*x) .* exp(-(x-0.5).^2 / 0.01);
gradu_Ex = @(x) exp(-100*(x - 1/2).^2).*sin(12*pi*x).*(200*x - 100) - 12*pi*exp(-100*(x - 1/2).^2).*cos(12*pi*x);
f = @(x) exp(-100*(x - 1/2).^2).*sin(12*pi*x).*(200*x - 100).^2 - 200*exp(-100*(x - 1/2).^2).*sin(12*pi*x) - ...
    144*pi^2*exp(-100*(x - 1/2).^2).*sin(12*pi*x) - 24*pi*exp(-100*(x - 1/2).^2).*cos(12*pi*x).*(200*x - 100);
k = @(x) 1 * x.^0;

fig = createFigure(W, H, 'enhanced',enhanced);
xPlot = linspace(0, 1, 1000);
plot(xPlot, uEx(xPlot), 'k')
xlabel('$x$', 'interpreter', 'laTeX')
legend({'$u(x)$'}, 'location', 'NE', 'interpreter', 'laTeX')
box on
axPos = get(gca, 'innerposition');
% export_fig(fig, '../../../Reports/DraftPDE_18/VERSION7/Figures/ErrEst1D_ExSol.eps', '-nocrop')
% close

p = 1;

NVec = 2.^[2:6];
hVec = 1 ./ NVec;

M = 10;

xFine = 0:1e-3:1;

for i = 1 : length(hVec)
    
    h = hVec(i);
    x = 0:h:1;
    N = length(x)-1;
    
    err = zeros(1, N);
    errProb = zeros(1, N);
    
    F = assembleRHS(f, x);
    A = assembleMatrix(k, x);
    u = [0; A \ F; 0];
    
    for kk = 1 : N
        err(kk) = sqrt(integral(@(xx) (uEx(xx) - interp1(x, u, xx)).^2, x(kk), x(kk+1)));
    end
    
    if i == 4
        fig = createFigure(W, H, 'enhanced', enhanced);
        axDP = axes;
        hold on
    end
    
    diff = 0;
    
    for jj = 1 : M
        
        X = x;
        for kk = 2 : length(x)-1
            hBar = min(x(kk)-x(kk-1), x(kk+1)-x(kk));
            X(kk) = x(kk) + hBar^p * (rand - 0.5);
        end
        
        FTilde = assembleRHS(f, X);
        ATilde = assembleMatrix(k, X);
        uTilde = [0; ATilde \ FTilde; 0];
        
        for kk = 1 : N
            errProb(kk) = errProb(kk) + 1/M * sqrt(integral(@(xx) (interp1(x, u, xx) - interp1(X, uTilde, xx)).^2, x(kk), x(kk+1)));
        end
        
        if i == 4 && jj == 1
            p1 = plot(axDP, X, uTilde, 'color', 0.7 * ones(3, 1));
        elseif i == 4
            plot(axDP, X, uTilde, 'color', 0.7 * ones(3, 1));
        end
        
        diff = diff + 1/(h*M) * abs(interp1(x, u, xFine) - interp1(X, uTilde, xFine));
        
    end
    
    if i == 4
        p2 = plot(axDP, x, u, 'k');
        legend(axDP, [p2 p1], {'$u_h$', '$U_h$'}, 'location', 'NE', 'interpreter', 'laTeX')
        box on
        xlabel(axDP, '$x$', 'interpreter', 'LaTeX')
        set(gca, 'innerposition', axPos)
%         export_fig(fig, '../../../Reports/DraftPDE_18/VERSION7/Figures/ErrEst1D_DetVsProb.eps', '-nocrop')
%         close
    end
    
    if i == 2
        fig = createFigure(W, H, 'enhanced', enhanced);
        semilogy((x(2:end) + x(1:end-1))/2, err, 'k.-');
        hold on
        semilogy((x(2:end) + x(1:end-1))/2, errProb, '.-', 'color', 0.7 * ones(3, 1));
        xlabel('$x$', 'interpreter', 'laTeX')
        legend({'true', 'est'}, 'location', 'NE', 'interpreter', 'laTeX')
        box on
        set(gca, 'innerposition', axPos)
%         export_fig(fig, '../../../Reports/DraftPDE_18/VERSION7/Figures/ErrEst1D_ErrCoarse.eps', '-nocrop')
%         close
    end
    
    if i == length(hVec)
        fig = createFigure(W, H, 'enhanced', enhanced);
        semilogy((x(2:end) + x(1:end-1))/2, err, 'k.-');
        hold on
        semilogy((x(2:end) + x(1:end-1))/2, errProb, '.-', 'color', 0.7 * ones(3, 1));
        xlabel('$x$', 'interpreter', 'laTeX')
        legend({'true', 'est'}, 'location', 'NE', 'interpreter', 'laTeX')
        box on
        set(gca, 'innerposition', axPos)
%         export_fig(fig, '../../../Reports/DraftPDE_18/VERSION7/Figures/ErrEst1D_ErrFine.eps', '-nocrop')
%         close
    end
    
end