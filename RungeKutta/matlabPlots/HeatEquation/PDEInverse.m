clc; clear; close all;

%% INVERSE PROBLEM PARABOLIC TEST

%% FINE-MESH MODEL

dxRef = 0.01;
T = 0.01;
xRef = 0 : dxRef : 1;
NRef = length(xRef) - 1;
hRef = 0.001;

ARef = -gallery('tridiag', NRef-1);
ARef = ARef / (dxRef^2);

% u0 = exp(-0.5 * (xRef - 0.5).^2 / (2 * 0.01));
u0 = (xRef < 0.9) .* (xRef > 0.7) + (xRef < 0.3) .* (xRef > 0.1);
u0 = u0';
plot(xRef, u0)

% Forward operator of the heat equation
forwMatrix = expm(ARef * T);

% Observation noise distribution
nEig = 1e4;
[V, D] = laplEig(xRef, nEig);
delta = 1e-1;
sqrtD = delta ./ diag(sqrt(-D));


% Generate noise and plot data
noise = KL(zeros(size(xRef')), V, sqrtD, randn(nEig, 1));
noise = noise(2:end-1);
yObs = forwMatrix * u0(2:end-1) + noise;

figure
plot(xRef(2:end-1), yObs)

%% INFERENCE SET-UP

% Compute EigenFunctions for inference
dx = 0.01;
x = 0 : dx : 1;
x = x';
N = length(x) - 1;
nKL = 100;
[V, D] = laplEig(x, nKL);
sqrtD = 1 ./ diag(sqrt(-D));
uMean = zeros(size(x));

% One sample from the prior
sample = KL(uMean, V, sqrtD, randn(nKL, 1));
plot(x, sample)

A = -gallery('tridiag', N-1);
A = A / (dx^2);
A = 1e-2 * A;

% Integration parameter
h = 0.01;

%% MCMC

% Data
nMCMC = 1e4;

% Do MCMC (CN)
posteriorCN = @(theta) likelihood_CN(theta(2:end-1), A, [0, T], h, x(2:end-1)', yObs', delta, x);
beta = 0.01;
proposalCN = @(theta) sqrt(1-beta^2) * theta + beta * KL(uMean, V, sqrtD, randn(nKL, 1));
[thetaAllCN, accRatioCN] = CNMetHas(sample, posteriorCN, nMCMC, proposalCN);
MCMCSample = thetaAllCN(:, ceil(nMCMC / 10) : end);

%% Summary plot

mcmcEst = mean(MCMCSample, 2);
mcmcVar = var(MCMCSample, 0, 2);

plot(xRef, u0);
hold on
plot(x, mcmcEst);

maxLimKappaPrc = prctile(MCMCSample', 97.5);
minLimKappaPrc = prctile(MCMCSample', 2.5);
confInt = [minLimKappaPrc, fliplr(maxLimKappaPrc)];
fill([x', fliplr(x')], confInt, 'b', ...
    'linestyle', '--', 'facealpha', 0.1);
legend('true', 'average', 'confidence 0.95', 'Location', 'best')

pU = figure;
close
% avgKappa = buildField_Andrea(mean(thetaAll, 2));

