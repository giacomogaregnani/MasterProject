clc; clear; close all;

%% INVERSE PROBLEM 2D TEST

%% CREATE FINE MESH

% Space discretization
hRef = 0.01;
[vertices, boundaries, elements] = initmesh('mesh_square','Jiggle','minimum','Hgrad',1.01,'Hmax', hRef);

%% FINE-MESH MODEL

param = [];
data   = read_DataFile('InverseData', 2, param);
data.diffusion = @(x, y, t, param)(1.3 + 0.5 * ((x - 0.25).^2 + (y - 0.25).^2 < 0.04) ...
                                       - 0.5 * ((x - 0.75).^2 + (y - 0.75).^2 < 0.04));
data.force = @(x,y,t,param)(8*pi^2*(sin(2*pi*x).*sin(2*pi*y)));
data.param = param;
data.bcDir = @(x,y,t,param)(0*x.*y);  

% Solve
[uRef, ~, meshRef]  = Elliptic_Solver(2, elements, vertices, boundaries, 'P1', data);

% Plot true solution
FEMplot2D(meshRef, uRef)
title('true solution')

% Plot true diffusion
diffNodes = data.diffusion(meshRef.nodes(1, :), meshRef.nodes(2, :), [], []);
FEMplot2D(meshRef, diffNodes)
title('true diffusion')

% Generate observations
[XX, YY] = meshgrid(linspace(0.1, 0.9, 6));
xObs = [XX(:)'; YY(:)'];

obsNoise = 1e-2;
uObs = evalFEM(meshRef, uRef, xObs);
uObs = uObs + obsNoise * randn(size(uObs));

%% CREATE COARSE MESH    

hCoarse = 0.2;
[vert, bound, el] = initmesh('mesh_square','Jiggle','minimum','Hgrad',1.01,'Hmax', hCoarse);

%% COMPUTE EIGENFUNCTIONS ON COARSE MESH

nKL = 40;
[uCoarse, ~, meshCoarse, ~]  = Elliptic_Solver(2, el, vert, bound, 'P1', data);
[V, D] = eigenFuncts(meshCoarse, nKL); 
sqrtD = diag(sqrt(diag(D)));

% One sample from the prior 
uMean = zeros(size(uCoarse));
sample = KL(uMean, V, sqrtD, randn(size(D, 1), 1));
FEMplot2D(meshCoarse, exp(sample));
title('sample from the prior')

%% MCMC

% prior choice
prior = @(x) -0.5 * x' * x;

% Initial guess
t0 = zeros(nKL, 1);

% Data
nMCMC = 5e4;

% Likelihood model
likelihood = @(theta) likelihood_FEM(theta, V, sqrtD, xObs, uObs, obsNoise, meshCoarse, data);
likelihoodProb = @(theta) likelihood_FEMProb(theta, V, sqrtD, xObs, uObs, obsNoise, meshCoarse, data);

% Posteriors
posterior = @(theta) prior(theta) + likelihood(theta);
posteriorProb = @(theta) prior(theta) + likelihoodProb(theta);

% Do MCMC (det)
sigma = min(meshCoarse.h);
[thetaAll, accRatio, S] = MetropolisHastings(t0, posterior, nMCMC, true, sigma);
thetaAllEffDet = thetaAll(:, ceil(nMCMC / 10) : end);

% Do MCMC (prob)
sigma = min(meshCoarse.h);
thetaAll = MetropolisHastings(t0, posteriorProb, nMCMC, true, sigma);
thetaAllEffProb = thetaAll(:, ceil(nMCMC / 10) : end);

% Results
MCMCEstDet = plotMCMCResultsPDE(thetaAllEffDet, uMean, V, sqrtD, meshCoarse);
MCMCEstProb = plotMCMCResultsPDE(thetaAllEffProb, uMean, V, sqrtD, meshCoarse);