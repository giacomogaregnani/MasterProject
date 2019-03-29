clc; clear; close all;
%%

thetaEx = [0.2, 0.2, 3.0];

plotParam = 3;

enhanced = 1;
fontsizeLAB = getLatexTextSize('normalsize', 'enhanced', enhanced);
fontsizeTICK = getLatexTextSize('small', 'enhanced', enhanced);

W = 4; H = 4;
figProb = createFigure(W, H, 'enhanced',enhanced);
mF = 0;

suff = {'1', '2', '3', '4', '5'};
varMC = [];
style = {'-', '-.', '--', '--.', ':'};
for j = 1 : length(suff)
    x = dlmread(['MCMC_ODE_PROB' suff{j} '.txt']);
    varMC = [varMC var(x(:, 3))];
    
    [f, xi] = ksdensity(exp(x(:, plotParam)), 'bandwidth', 0.03);

    plot(xi, f, style{j}, 'markersize', 3.5)
    hold on
    
    mF = max(mF, max(f));
    
end
xlim([2.5, 4])
xLim = get(gca, 'xLim');
xlabel('$e^{\vartheta_3}$', 'interpreter', 'LaTex')
plot(figProb.CurrentAxes, [thetaEx(plotParam), thetaEx(plotParam)], [0, mF], 'k--')

set(gca, 'yTick', []);
set(gca, 'fontsize', fontsizeTICK);
set(get(gca, 'xlabel'), 'fontsize', fontsizeLAB);
set(get(gca, 'ylabel'), 'fontsize', fontsizeLAB);
box on
leg = {'h = h_0', 'h = h_1', 'h = h_2', 'h = h_3', 'h = h_4'};
% legend(leg, 'location', 'best')

fig = createFigure(W, H, 'enhanced',enhanced);
mF = 0;
for j = 1 : length(suff)
    x = dlmread(['MCMC_ODE_DET' suff{j} '.txt']);
    
    [f, xi] = ksdensity(exp(x(:, plotParam)), 'bandwidth', 0.03);
    
    plot(xi, f, style{j}, 'markersize', 3.5)
    hold on
    
    mF = max(mF, max(f));
    
end
xlabel('$e^{\vartheta_3}$', 'interpreter', 'LaTex')

set(gca, 'xLim', xLim)
set(gca, 'yTick', []);
set(gca, 'fontsize', fontsizeTICK);
set(get(gca, 'xlabel'), 'fontsize', fontsizeLAB);
set(get(gca, 'ylabel'), 'fontsize', fontsizeLAB);
yLim = get(gca, 'yLim');
box on
% legend(leg, 'location', 'best')
plot(fig.CurrentAxes, [thetaEx(plotParam), thetaEx(plotParam)], [0, mF], 'k--')
print -depsc2 ../../../Reports/PaperRTSRK_18/VERSION11/PostDet.eps
% close(fig)

set(figProb.CurrentAxes, 'yLim', yLim)

print -depsc2 ../../../Reports/PaperRTSRK_18/VERSION11/PostProb.eps

%% Test two dimensional param

% thetaProb = dlmread('MCMC_ODE_PROB1.txt');
% thetaDet = dlmread('MCMC_ODE1.txt');
% 
% i = 1; j = 3;
% 
% nPoints = 50;
% xEval = linspace(min(thetaProb(:, 1)), max(thetaProb(:, 1)), nPoints);
% yEval = linspace(min(thetaProb(:, 3)), max(thetaProb(:, 3)), nPoints);  
% [XX, YY] = meshgrid(xEval, yEval);
% 
% dens = akde([thetaProb(:, 1), thetaProb(:, 3)], [XX(:), YY(:)]);
% dens = reshape(dens, nPoints, nPoints);
% 
% xEvalDet = linspace(min(thetaDet(:, 1)), max(thetaDet(:, 1)), nPoints);
% yEvalDet = linspace(min(thetaDet(:, 3)), max(thetaDet(:, 3)), nPoints);
% [XXDet, YYDet] = meshgrid(xEvalDet, yEvalDet);
% densDet = akde([thetaDet(:, 1), thetaDet(:, 3)], [XX(:), YY(:)]);
% densDet = reshape(densDet, nPoints, nPoints);
% % 
% %% Plots
% 
% figure
% contour(XX, YY, dens, 20)
% hold on
% contour(XXDet, YYDet, densDet, 20)
% plot(thetaEx(i), thetaEx(j), 'or')