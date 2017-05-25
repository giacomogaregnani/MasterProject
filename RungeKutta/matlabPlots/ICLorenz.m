clear; clc; 

addpath('ResultsICLorenz');

info = dlmread('infoFile.txt');
T = info(1);
h = info(2);
timeVec = 0 : h : T;

% Only final value, plot an histogram
solution = dlmread('ICLorenz.txt');
phiSolution = sum(solution.^2, 2);
histogram(phiSolution, 'normalization', 'pdf', 'facealpha' , .4)

% kernel density estimation

evalPoints = [solution(:, 1), solution(:, 2), solution(:, 3)];
V = normalKernelDensity(solution, evalPoints);

%% Plot
close all;
nPointsForPlot = 10000;

scatter3(solution(1:nPointsForPlot, 1), solution(1:nPointsForPlot, 2), solution(1:nPointsForPlot, 3), 10.0, V(1:nPointsForPlot))
colorbar
colormap('parula')

figure
scatter(solution(1:nPointsForPlot, 1), solution(1:nPointsForPlot, 3), 10.0, V(1:nPointsForPlot))
colormap('parula')
