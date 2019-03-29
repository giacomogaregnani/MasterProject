clc; clear; close all
%% Model

uEx = [];

f = @(x) sin(2 * pi * x);
field = @(x, param) 1  + exp(param(1)) * (x>=0.2) .* (x<=0.4) + exp(param(2)) * (x>=0.6) .* (x<=0.8);

trueParam = [0.3; 1];

% Generate observations

sigma = 1e-4;
xObs = linspace(0.1, 0.9, 20)';

if isempty(uEx)
    
    disp('Computing observations...');
    
    xEx = linspace(0, 1, 1e3+1);
    F = assembleRHS(f, xEx);
    A = assembleMatrix(@(x) field(x, trueParam), xEx);
    
    uExVec = [0; A \ F; 0];
    uEx = @(xx) interp1(xEx, uExVec', xx);
    
    observations = uEx(xObs) + sigma * randn(size(xObs));
    
    disp('Computed observations');
    
end

%% Inverse problem
NVec = [10, 20, 40, 80];

sampleDet = cell(length(NVec), 1);
sampleProb = sampleDet;

nChains = 50;
nMCMC = 2000;

prior = @(param) -0.5 * (param') * param;
proposal = @(param) gaussianProposal(param, 500 * sigma);

for i = 1 : length(NVec)
    x = linspace(0, 1, NVec(i));

    likelihoodDet  = @(param) likelihoodFEM(param, field, f, x, xObs, observations, sigma);
    posteriorDet = @(param) likelihoodDet(param) + prior(param);
    initGuess = trueParam;
    
    sampleDet{i} = metropolisHastings(posteriorDet, proposal, initGuess, nMCMC, false);
    sampleProb{i} = metropolisHastingsProbFEM(proposal, prior, x, field, f, xObs, observations, sigma, initGuess, nMCMC, nChains);
end

save('FinDimInverseProblem2')

%% Plot results

% One-dimensional plots
col = parula(2*length(NVec));
col = col(1:2:end, :);
for j = 1:2
    figure
    hold on
    for i = 1 : length(NVec)
        [dens, xi] = ksdensity(sampleDet{i}(j, :));
        plot(xi, dens, 'color', col(i, :));
        [dens, xi] = ksdensity(sampleProb{i}(j, :));
        plot(xi, dens, '--','color', col(i, :));
    end
    plot([trueParam(j), trueParam(j)], [0, 3*max(dens)], 'k--')
end

% Two-dimensional plots
figure
hold on
for i = 1 : length(NVec)
    plotTwoDimDens(sampleProb{i}', 50, col(i, :))
end
plot(trueParam(1), trueParam(2), 'xk')
xLim = get(gca, 'xlim'); yLim = get(gca, 'ylim');

figure
hold on
for i = 1 : length(NVec)
    plotTwoDimDens(sampleDet{i}', 50, col(i, :))
end 
plot(trueParam(1), trueParam(2), 'xk')
set(gca, 'xlim', xLim); set(gca, 'ylim', yLim);

    
%% Solution on true param and on mean param

% F = assembleRHS(f, x);
% A = assembleMatrix(@(xx) field(xx, trueParam), x);
% uVecEx = [0; A \ F; 0];
%
% meanParam = mean(sample,2);
% covParam = cov(sample');
%
% A = assembleMatrix(@(xx) field(xx, meanParam), x);
% uVecMean = [0; A \ F; 0];
%
% figure
% plot(xEx, uEx(xEx))
% hold on
% plot(x, uVecEx)
% plot(x, uVecMean)
% legend('exact', 'ex param', 'mean param')


