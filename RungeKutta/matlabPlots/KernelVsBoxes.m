%% Compare results of Kernel density and boxes
clc; clear; close all;

addpath('ResultsLorenzBoxes')
addpath('ResultsStep')

%% Read boxes file, compute density

% Parameters boxes
origin = [-30, -30, 0];
bigRadius = [60, 60, 60];

[densBoxes, ~, xBoxes] = computeDensity('testCubesD18H05', 1e-12);


%% Read step file, compute density on box edges

solutionStep = dlmread('SolutionStep.txt');
densStep = akde(solutionStep, xBoxes(:, 1 : 3) + xBoxes(:, 4 : 6)/2);

%% Plot results

transp = true;

plotBoxCollectionColor(xBoxes, densBoxes, transp);
colorLim = get(gca, 'CLim');
colorbar
title('Subdivision')

plotBoxCollectionColor(xBoxes, densStep, transp);
set(gca, 'CLim', colorLim);
colorbar
title('Random time-stepping')

%% Compute error

% volume of each element = product of the three radii
vol = prod(xBoxes(1, 4 : 6));

intL2Step = vol * sum(densStep.^2);
intL2Boxes = vol * sum(densBoxes.^2);
errL2Norm = vol * sum((densStep - densBoxes).^2);