clc; clear; close all;
addpath('resultsOxPerox');
%%
nMC = 50;
x = dlmread('OxPeroxStep.txt');
N = size(x, 1) / nMC;

%%
% enhanced = 1;
% fontsizeLAB = getLatexTextSize('normalsize', 'enhanced', enhanced);
% fontsizeTICK = getLatexTextSize('small', 'enhanced', enhanced);

%%

W = 12; H = 3;
fig = createFigure(W, H, 'enhanced',1);
for i = 1 : 20
    idx = [(i-1)*N + 1 : i*N];
    semilogy(x(idx, 1), x(idx, 4), 'color', [0.7, 0.7, 0.7]);
    hold on
    xlabel('$t$', 'interpreter', 'laTeX')
end
box on
xlim([0, 200])
ylim([1e-5, 1])
set(gca, 'yTick', [1e-5, 1])
ylabel('[X]', 'interpreter', 'latex')

print -depsc2 ../../../Reports/PaperRTSRK_18/VERSION11/OxPerox.eps