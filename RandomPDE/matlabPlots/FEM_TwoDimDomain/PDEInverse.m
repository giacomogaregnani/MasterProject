clc; clear; close all;

%% INVERSE PROBLEM 2D TEST

%% CREATE FINE MESH

% Space discretization
hRef = 0.01;
[vertices, boundaries, elements] = poimesh('squareg', 500, 500);
vertices = 0.5 * vertices + 0.5;
% [vertices, boundaries, elements] = initmesh('mesh_square','Jiggle','minimum','Hgrad',1.01,'Hmax', hRef);

%% FINE-MESH MODEL

param = [];
data   = read_DataFile('InverseData', 2, param);
data.diffusion = @(x, y, t, param)(1.3 + 0.5 * ((x - 0.25).^2 + (y - 0.25).^2 < 0.025) ...
    - 0.5 * ((x - 0.75).^2 + (y - 0.75).^2 < 0.025));
% data.diffusion = @(x, y, t, param) 1 + x.*y;
data.force = @(x,y,t,param)(0*x.*y);
data.param = param;
data.bcDir = @(x,y,t,param)(0*x.*y);
data.bcNeu = @(x,y,t,param)((x - 1).^2);

% Solve
[uRef, ~, meshRef]  = Elliptic_Solver(2, elements, vertices, boundaries, 'P1', data);

% Plot true solution
FEMplot2D(meshRef, uRef)
title('true solution')

% Generate observations
nObsPerSide = 15; obsVec = linspace(0.1, 0.9, nObsPerSide);
xObs = [obsVec, obsVec, ones(1, nObsPerSide);
        zeros(1, nObsPerSide), ones(1, nObsPerSide), obsVec];
% [XX, YY] = meshgrid(linspace(0.1, 0.9, 9));
% xObs = [XX(:)'; YY(:)']; 


obsNoise = 1e-3;
uObs = evalFEM(meshRef, uRef, xObs);
uObs = uObs + obsNoise * randn(size(uObs));

% Plot observations
hold on
plot3(xObs(1, :), xObs(2, :), uObs, 'ro')

% Plot true diffusion
diffNodes = data.diffusion(meshRef.nodes(1, :), meshRef.nodes(2, :), [], []);
% FEMplot2D(meshRef, diffNodes)
title('true diffusion')


%% CREATE COARSE MESH

for hCoarse = 0.1
    
%     [vert, bound, el] = initmesh('mesh_square','Jiggle','minimum','Hgrad',1.01,'Hmax', hCoarse);
    [vert, bound, el] = poimesh('squareg', 40, 40);
    vert = 0.5 * vert + 0.5;

    %% COMPUTE EIGENFUNCTIONS ON COARSE MESH
    
    [uCoarse, ~, meshCoarse, ~]  = Elliptic_Solver(2, el, vert, bound, 'P1', data);
    nKL = size(vert, 2);
    [V, D] = eigenFuncts(meshCoarse, nKL);
    sqrtD = diag(sqrt(diag(D)));
    
    % One sample from the prior
    uMean = zeros(size(uCoarse));
    sample = KL(uMean, V, sqrtD, randn(size(D, 1), 1));
    FEMplot2D(meshCoarse, exp(sample));
    title('sample from the prior')
    
    %% MCMC
    
%     prior choice
    prior = @(x) -0.5 * x' * x;
    
%     Initial guess
    t0 = zeros(nKL, 1);
    
%     Data
    nMCMC = 15e3;
    
%     Likelihood model
    likelihood = @(theta) likelihood_FEM(theta, uMean, V, sqrtD, xObs, uObs, obsNoise, meshCoarse, data);
%     likelihoodProb = @(theta) likelihood_FEMProb(theta, V, sqrtD, xObs, uObs, obsNoise, meshCoarse, data);
    
%     Posteriors
    posterior = @(theta) prior(theta) + likelihood(theta);
%     posteriorProb = @(theta) prior(theta) + likelihoodProb(theta);
    
%     Do MCMC (det)
    sigma = hCoarse / 4;
    [thetaAll, accRatio, S] = MetropolisHastings(t0, posterior, nMCMC, true, sigma);
    thetaAllEffDet = thetaAll(:, ceil(nMCMC / 10) : end);
    
%     Do MCMC (prob)
%     sigma = hCoarse / 4;
%     thetaAll = MetropolisHastings(t0, posteriorProb, nMCMC, true, sigma);
%     thetaAllEffProb = thetaAll(:, ceil(nMCMC / 10) : end);
    
%     Do MCMC (CN)
%     posteriorCN = @(theta) likelihood_FEM(theta, uMean, V, sqrtD, xObs, uObs, obsNoise, meshCoarse, data);
%     beta = 0.02;
%     proposalCN = @(theta) sqrt(1-beta^2) * theta + beta * randn(nKL, 1);
%     [thetaAllCN, accRatioCN] = CNMetHas(randn(nKL, 1), posteriorCN, nMCMC, proposalCN);
%     thetaAllEffCN = thetaAllCN(:, ceil(nMCMC / 10) : end);
    
%     Results
    MCMCEstDet = plotMCMCResultsPDE(thetaAllEffDet, uMean, V, sqrtD, meshCoarse);
%     MCMCEstProb = plotMCMCResultsPDE(thetaAllEffProb, uMean, V, sqrtD, meshCoarse);
%     MCMCEstCN = plotMCMCResultsPDE(thetaAllEffCN, uMean, V, sqrtD, meshCoarse);
%     close all
%     save(['results_h_', num2str(hCoarse * 1000), '2.mat']);
    
end