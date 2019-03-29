clc; clear; close all
%%

addpath('resultsErrorEstimators')
enhanced = 1;
fontsizeLAB = getLatexTextSize('normalsize', 'enhanced', enhanced);
fontsizeTICK = getLatexTextSize('small', 'enhanced', enhanced);
W = 12; H = 4;

problem = 3;

switch problem
    case 1
        resultsProb = dlmread('fixedhFitzNag.txt');
        resultsDet = dlmread('fixedhFitzNagdet.txt');
        d = 2;
    case 2
        resultsProb = dlmread('fixedhHenHeil.txt');
        resultsDet = dlmread('fixedhHenHeildet.txt');
        d = 4;
    case 3
        resultsProb = dlmread('fixedhLorenz.txt');
        resultsDet = dlmread('fixedhLorenzdet.txt');
        d = 3;
    case 4
        resultsProb = dlmread('fixedhVdPol.txt');
        resultsDet = dlmread('fixedhVdPoldet.txt');
        d = 2;
end

t = resultsProb(:, 1);
h = t(2) - t(1);

c = colormap(parula(3));
close


fig = createFigure(W, H, 'enhanced',enhanced);
% subsemilogy(3, 1, 1)
semilogy(t(1:4000), resultsDet(1:4000, 3), 'k')
hold on
semilogy(t(1:4000), resultsDet(1:4000, 2), 'b')
semilogy(t(1:4000), resultsProb(1:4000, 2), 'r')
% plot(t, resultsProb(:, 2), 'r')
% plot(t, resultsProb(:, 5), 'c')
legend({'true error', 'embedded estimator', 'probabilistic estimator'}, 'interpreter', 'latex', 'location', 'SE')
axpos = get(gca, 'innerposition');
xlabel('$t$', 'interpreter', 'latex')
ylim([1e-4 1e2])
% print -depsc2 ../../../Reports/PaperRTSRK_18/VERSION11/aPosteriori.eps %-nocrop


% fig = createFigure(W, H, 'enhanced',enhanced);
% % subsemilogy(3, 1, 2)
% semilogy(t, resultsProb(:, 2), 'r')
% hold on
% semilogy(t, resultsProb(:, 3), 'b')
% % semilogy(t, resultsProb(:, 4), 'k')
% semilogy(t, resultsDet(:, 3), 'k')
% ylim([1e-4, 1e2])
% legend({'$tr(Var(Y))^{1/2}$', 'MSE', 'error'}, 'location', 'SE', 'interpreter', 'latex')
% set(gca, 'innerposition', axpos)
% % export_fig ../../../Reports/ReportPOST_18/VarianceMSE2.eps -nocrop
% 
% fig = createFigure(W, H, 'enhanced',enhanced);
% % subsemilogy(3, 1, 3)
% plot(t, resultsProb(:, 2) ./ resultsProb(:, 3), 'b')
% hold on
% plot(t, resultsProb(:, 2) ./ resultsDet(:, 3), 'k')
% legend({'Var / MSE', 'Var / error'}, 'interpreter', 'latex', 'location', 'NW')
% ylim([0, 5])
% set(gca, 'innerposition', axpos)
% % export_fig ../../../Reports/ReportPOST_18/VarianceErrorRatioVdPol2.eps -nocrop
% 
% % figure
% % semilogy(t, resultsDet(:, 4))
% % hold on
% % semilogy(t, resultsDet(:, end - d + 1))