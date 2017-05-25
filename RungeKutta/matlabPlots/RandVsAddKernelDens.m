clear; clc; close all;

%% Load numerical solutions and info
addpath('ResultsAdd', 'ResultsStep');

solutionAddFile = 'SolutionAdd.txt';
solutionStepFile = 'SolutionStep.txt';
solutionAdd = dlmread(solutionAddFile);
solutionStep = dlmread(solutionStepFile);

info = dlmread('infoFile.txt');
nMC = info(1);
T = info(2);
h = info(3);
fullWrite = info(4);
timeVec = 0 : h : T;

%% Kernel density estimation
    
% VAdd = normalKernelDensity(solutionAdd, solutionAdd);
VAdd = akde(solutionAdd, solutionAdd);
% VStep = normalKernelDensity(solutionStep, solutionStep);
VStep = akde(solutionStep, solutionStep);

%% Plot
close all

lines = dlmread('SolutionTimePlot.txt');

scatter3(solutionStep(:, 1), solutionStep(:, 2), ...
    solutionStep(:, 3), 8.0, VStep(:), 'filled')
colorbar
colormap('parula')
view([57, 15])
colorLim = get(gca, 'CLim');

hold on
% plot3(lines(:, 1), lines(:, 2), lines(:, 3), 'color', [.3, .3, .3]);

% 
% figure
% scatter(solutionStep(:, 1), solutionStep(:, 3), 10.0, VStep(:))
% colormap('parula')
% colorbar

figure
scatter3(solutionAdd(:, 1), solutionAdd(:, 2), ...
    solutionAdd(:, 3), 8.0, VAdd(:), 'filled')
colorbar
colormap('parula')
view([57, 15])
set(gca, 'CLim', colorLim);

hold on
% plot3(lines(:, 1), lines(:, 2), lines(:, 3), 'color', [.3, .3, .3]);

% figure
% scatter(solutionAdd(:, 1), solutionAdd(:, 3), 10.0, VAdd(:))
% colormap('parula')
% colorbar