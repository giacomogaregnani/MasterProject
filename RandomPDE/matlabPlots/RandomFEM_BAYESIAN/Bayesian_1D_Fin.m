clc; clear; close all
%%

N = 20;
for i = 1 : 3
    N(i+1) = N(i) * 2;
end
trueParam = [0.3, -0.3];

col = [0.00, 0.45, 0.74;  % blu
       0.49, 0.18, 0.56;  % viola
       0.85, 0.33, 0.10;  % rosso
       0.00, 0.50, 0.00]; % verde

W = 6.7; H = 6.7;
fontsizeLAB = getLatexTextSize('normalsize', 'enhanced', 1);
fontsizeTICK = getLatexTextSize('small', 'enhanced', 1);
   
fig = createFigure(W, H, 'enhanced', 1);
hold on
for i = 1 : length(N)
    sampleProb = dlmread(['finDim'  num2str(N(i)) 'coeffsProb.txt']); 
    plotTwoDimDens(sampleProb, 40, col(i, :))
end
plot(trueParam(1), trueParam(2), 'xk', 'linewidth', 1)
xLim = [ 0.2, 0.8];
yLim = [-0.4, 0  ];
set(gca, 'xlim', xLim); set(gca, 'ylim', yLim);
% xLim = get(gca, 'xlim'); yLim = get(gca, 'ylim');
legText = cellfun(@(c)['$N = ' c '$'],strsplit(num2str(N)),'uni',false);
legend(legText, 'interpreter', 'latex', 'position', [0.5417 0.6873 0.3545 0.2275])
xlabel('$\kappa_1$', 'interpreter', 'latex')
ylabel('$\kappa_2$', 'interpreter', 'latex')
box on
axPos = get(gca, 'innerposition');
% export_fig(fig, '../../../Reports/DraftPDE_18/VERSION7/Figures/Bayes_Fin_Prob.eps', '-nocrop')

fig = createFigure(W, H, 'enhanced',1);
hold on
for i = 1 : length(N)
    sampleDet = dlmread(['finDim'  num2str(N(i)) 'coeffs.txt']); 
    plotTwoDimDens(sampleDet, 40, col(i, :))
end
plot(trueParam(1), trueParam(2), 'xk', 'linewidth', 1)
legend(legText, 'interpreter', 'latex', 'position', [0.5417 0.6873 0.3545 0.2275])
set(gca, 'xlim', xLim); set(gca, 'ylim', yLim);
xlabel('$\kappa_1$', 'interpreter', 'latex')
ylabel('$\kappa_2$', 'interpreter', 'latex')
box on
set(gca, 'innerposition', axPos)
% export_fig(fig, '../../../Reports/DraftPDE_18/VERSION7/Figures/Bayes_Fin_Det.eps', '-nocrop')