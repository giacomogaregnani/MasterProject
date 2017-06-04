clc; close all; clear

%% Load data

solDet = dlmread('LorenzDet.txt');
solProb = dlmread('LorenzProb.txt');

% get Parameters
info = dlmread('LorenzProbInfo.txt');
h = info(1); T = info(2); nMC = info(3);

%% Plot trajectories

% Time vector
t = 0 : h : T;

lightgrey = [0.6, 0.6, 0.6];

enhanced = 1;
W = 10;
H = 3;

fontsizelb = 'normalsize';
fontsizetk = 'small';

fig = createFigure(W, H, 'enhanced',enhanced);
fontsizeLAB = getLatexTextSize(fontsizelb, 'enhanced', enhanced);
fontsizeTICK = getLatexTextSize(fontsizetk, 'enhanced', enhanced);

for i = 1 : nMC
    index = (i-1)*length(t)+1;
    plot(t, solProb(index:index+length(t)-1, 1), 'color', lightgrey)
    hold on
end

plot(t, solDet(:, 1), 'k', 'LineWidth', 1)

print -depsc2 ../../../Reports/Leysin2017/Lorenz.eps
