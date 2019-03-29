clc; clear; close all;

%%
nMC = 100;
xDet = dlmread('testVardet.txt');
x = dlmread('testVar.txt');
N = size(x, 1) / nMC;

%%
% enhanced = 1;
% fontsizeLAB = getLatexTextSize('normalsize', 'enhanced', enhanced);
% fontsizeTICK = getLatexTextSize('small', 'enhanced', enhanced);

%%

% W = 15; H = 5;
% fig = createFigure(W, H, 'enhanced',enhanced);
figure
hold on
for i = 1 : 20
    idx = [(i-1)*N + 1 : i*N];
    plot(x(idx, 1), x(idx, 2), 'color', [0.7, 0.7, 0.7]);
    xlabel('$t$', 'interpreter', 'laTeX')
end
plot(xDet(:, 1), xDet(:, 2), 'k', 'linewidth', 1)
box on

% print -depsc2 ../../../Reports/RandTimeStep/VERSION7/Lorenz.eps


%%

% times = [10 20 30];
% for i = times
%     
%     idx = find(x(:, 1) == i);
%     val = x(idx, 2);
%     
%     W = 4; H = 4;
%     fig = createFigure(W, H, 'enhanced',enhanced);
%     
%     [f, xi] = ksdensity(val);
%     plot(xi, f, 'k');
%     set(gca, 'yTick', [])
%     
%     if i == times(3)
%         ylim = get(gca, 'yLim');
%     end
%     if i > times(3)
%         set(gca, 'yLim', ylim);
%     end
%     box on
%     
%     set(gca, 'fontsize', fontsizeTICK);
%     set(get(gca, 'xlabel'), 'fontsize', fontsizeLAB);
%     set(get(gca, 'ylabel'), 'fontsize', fontsizeLAB);
%    
%     xlim([-25, 25])
%     view([90 -90])
%     
%     yTextPlot = get(gca, 'yLim');
%     text(20, yTextPlot(2) / 2, ['$t = $' num2str(i)], 'interpreter', 'laTeX')
%     
%     filename = [ '../../../Reports/RandTimeStep/VERSION7/Lorenz' num2str(i), '.eps'];
%     
% %     print('-deps', filename)
%     
% end


