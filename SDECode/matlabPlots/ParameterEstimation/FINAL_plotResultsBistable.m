clc; clear; close all
addpath('resultsBistable_GOOD')
%%

homParam = dlmread('ResultsHom.txt');
aHom = homParam(1:2);

W = 4.5; H = 4.5;
enhanced = 1;
fontsizeLAB = getLatexTextSize('small', 'enhanced', 1);
fontsizeTICK = getLatexTextSize('scriptsize', 'enhanced', 1);

for T = [100, 200, 400]
    
    aMLE = dlmread(['ResultsFilterT_' num2str(T) 'MLE.txt']);
    aMean = dlmread(['ResultsFilterT_' num2str(T) 'Mean.txt']);
    aCov = dlmread(['ResultsFilterT_' num2str(T) 'Cova.txt']);
    
    nPoints = 200;
    [X1, X2] = meshgrid(linspace(aMean(1)-3*sqrt(aCov(1,1)), aMean(1)+3*sqrt(aCov(1,1)), nPoints), ...
        linspace(aMean(2)-3*sqrt(aCov(2,2)), aMean(2)+3*sqrt(aCov(2,2)), nPoints));
    X = [X1(:), X2(:)];
    p = mvnpdf(X, aMean', aCov);
    
    fig = createFigure(W, H, 'enhanced', enhanced);
    contour(X1, X2, reshape(p, nPoints, nPoints), 'color', 0.7 * ones(1, 3));
    hold on
    plot(aMLE(1), aMLE(2), 'ko', 'markersize', 5)
    plot(aHom(1), aHom(2), 'kx', 'markersize', 5)
    legend({'$\mu$', '$\widetilde A^{\varepsilon}_k$', '$A$'}, 'location', 'nw', 'interpreter', 'latex')
    
    xlim([0.1, 0.52])
    ylim([0.25, 1.2])
    
    xlabel('$A_1$', 'interpreter', 'latex')
    ylabel('$A_2$', 'interpreter', 'latex')
        
    title(['$T = $ ', num2str(T)], 'interpreter', 'latex');
    
    set(gca, 'fontsize', fontsizeTICK);
    set(get(gca, 'xlabel'), 'fontsize', fontsizeLAB);
    set(get(gca, 'ylabel'), 'fontsize', fontsizeLAB);
    set(get(gca, 'title'), 'fontsize', fontsizeLAB);
    set(get(gca, 'legend'), 'fontsize', fontsizeLAB);
    
    placeLegendExtreme(fig, 3, 'NW')
    niceBox(fig)
    
    export_fig(fig, ['~/Desktop/PaperSDE/Figures/Bayes' '_T' num2str(T) '.png'], '-nocrop', '-painters', '-m5')
end