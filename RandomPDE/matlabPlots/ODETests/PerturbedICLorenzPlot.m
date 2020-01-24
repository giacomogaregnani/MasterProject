clc; clear; close all;
addpath('resultsPertICLorenz');
%%
nMC = 20;
x = dlmread('pertICLorenz1e-5.txt');
N = size(x, 1) / nMC;

%%
% enhanced = 1;
% fontsizeLAB = getLatexTextSize('normalsize', 'enhanced', enhanced);
% fontsizeTICK = getLatexTextSize('small', 'enhanced', enhanced);

%%

W = 12; H = 3;
fig = createFigure(W, H, 'enhanced',1);
hold on
for i = 1 : 20
    idx = (i-1)*N+1 : 10 : i*N;
    plot(x(idx, 1), x(idx, 2), 'color', [0.7, 0.7, 0.7]);
    xlabel('$t$', 'interpreter', 'laTeX')
    ylabel('$y_1(t)$', 'interpreter', 'laTeX')
end
box on

% print -depsc2 ../../../Reports/PaperRTSRK_18/VERSION17/Paper/LorenzTest3.eps


