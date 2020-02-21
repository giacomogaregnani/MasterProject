clc; clear; close all
addpath('resultsBayesianOne')
%%

homParam = dlmread('ResultsHom.txt');
aHom = homParam(1);

W = 4; H = 4;
enhanced = 4;
fontsizeLAB = getLatexTextSize('footnotesize', 'enhanced', 1);
fontsizeTICK = getLatexTextSize('scriptsize', 'enhanced', 1);

errMeanMLE = zeros(5, 1);
times = [100, 200, 300, 400, 500];

% for i = 1 : length(times)

T = 2000; % times(i);

aMLE = dlmread(['ResultsFilterT2_' num2str(T) 'MLE.txt']);
aMean = dlmread(['ResultsFilterT2_' num2str(T) 'Mean.txt']);
aVar = dlmread(['ResultsFilterT2_' num2str(T) 'Cova.txt']);

%     nPoints = 1000;
%     X = linspace(aMean-4*sqrt(aVar), aMean+4*sqrt(aVar));
%     p = normpdf(X, aMean, sqrt(aVar));
%     maxp = max(p);
%
%     fig = createFigure(W, H, 'enhanced', enhanced);
%     plot(X, p, 'k');
%     hold on
%     plot([aMLE, aMLE], [0, maxp], 'k--')
%     plot([aHom, aHom], [0, maxp], 'k-.')
%     legend({'$\mu$', '$\widehat A$', '$A$'}, 'location', 'nw', 'interpreter', 'latex')
%
%     xlim([0, 1])
%     ylim([0,12])
%
% %     xlabel('$A_1$', 'interpreter', 'latex')
% %     ylabel('$A_2$', 'interpreter', 'latex')
%
%     title(['$t = $ ', num2str(T)], 'interpreter', 'latex');
%
%     set(gca, 'fontsize', fontsizeTICK);
%     set(get(gca, 'xlabel'), 'fontsize', fontsizeLAB);
%     set(get(gca, 'ylabel'), 'fontsize', fontsizeLAB);
%     set(get(gca, 'title'), 'fontsize', fontsizeLAB);
%
%     export_fig(fig, ['~/Desktop/PaperSDE/Figures/Bayes' '_T' num2str(T) '.png'], '-nocrop', '-painters', '-m5')
% end

%% PLOTS
h = 0.05^3;
timeVec = 0 : h : T;
plot(timeVec, [0; aMean], 'k')
hold on
plot([0, T], [aMLE, aMLE], 'k--')
figure
semilogy(timeVec, [aMLE; abs(aMLE - aMean)]);
hold on
semilogy(timeVec, 1./timeVec);
ylim([0, 1])

figure
semilogy(timeVec, [1; aVar], 'k')
hold on
semilogy(timeVec, 1./timeVec);
ylim([0, 1])

%% ANIMATION
figure
maxp = 20;
for i = 1 : 8000 : length(timeVec)-1
    X = linspace(aMean(i)-4*sqrt(aVar(i)), aMean(i)+4*sqrt(aVar(i)));
    p = normpdf(X, aMean(i), sqrt(aVar(i)));
    
    plot([aMLE, aMLE], [0, maxp], 'k--')
    hold on
    plot([aHom, aHom], [0, maxp], 'k-.')
    plot(X, p, 'k');
    xlim([0, 1])
    ylim([0, 20])
    title(['t = ', num2str(timeVec(i))])
    pause(0.001)
    hold off
end




