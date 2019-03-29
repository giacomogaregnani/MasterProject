clc; clear; close all
%%

run('../initialize.m') % Insert path to the script 'initialize.m'

W = 10; H = 5; % Insert width and heigth of the figure 

fontsizeLAB = getLatexTextSize('normalsize', 'enhanced', 1);
fontsizeTICK = getLatexTextSize('small', 'enhanced', 1);

t = 0 : 0.01 : 1;
y = sin(2*pi*t);
   
fig = createFigure(W, H);
hold on
plot(t, y, 'k')
% plot(t, z, 'k--')
xlabel('$t$', 'interpreter', 'latex', 'fontsize', fontsizeLAB)
ylabel('$f(t, \varepsilon)$', 'interpreter', 'latex', 'fontsize', fontsizeLAB)
title({'Some sinusoidals'}, 'interpreter', 'latex', 'fontsize', fontsizeLAB)
box on
export_fig(fig, 'test.eps', '-nocrop') % This figure should be included with \includegraphics[]{test} in your .tex file