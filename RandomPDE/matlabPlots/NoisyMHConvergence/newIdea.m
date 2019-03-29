clc; clear; close all
%% Generate observations

thetaEx = 1;
f = @(x) sin(2 * pi * x)';
xRef = linspace(0, 1, 10000);
uEx = @(x) 1 / (2 * pi)^2 * sin(2 * pi * x);

% Generate observations
xObs = 0.4;
obsNoise = 1e-4;
uObs = 1 / exp(thetaEx) * uEx(xObs) + obsNoise * randn(1);

%% Inference

nMCMC = 2e4;
sigmaProp = 0.1;
sigmaPropProb = sigmaProp * 20;
proposal = @(w) GaussianProposal(w, sigmaProp);
priorMean = 0;
priorVariance = 1;
prior = @(w) GaussianPrior(w, priorMean, priorVariance);

% Exact inference
likelihoodEx = @(w) ExactLikelihood(w, xObs, uObs, obsNoise);
posteriorEx = @(w) likelihoodEx(w) + prior(w);
theta = MetropolisHastings(priorMean, posteriorEx, proposal, nMCMC);
theta = theta(floor(nMCMC / 10) : end);
exEstimate = mean(theta);

% FEM-based inference
NVec = [1, 2, 4, 8, 16];
estDist = [];
estDistProb = [];
N = 10;
h = 1 / N;
x = linspace(0, 1, N+1);
A = 1 / h * gallery('tridiag', N-1);
F = h * f(x(2:end-1));

for nMC = NVec
    
%     likelihoodFEM = @(w) FEMLikelihood(w, xObs, uObs, obsNoise, A, F, x);
%     posteriorFEM = @(w) prior(w) + likelihoodFEM(w);
%     thetaFEM = MetropolisHastings(priorMean, posteriorFEM, proposal, nMCMC);
%     thetaFEM = thetaFEM(floor(nMCMC / 10) : end);
%     estDist = [estDist abs(mean(thetaFEM) - exEstimate)];
  
    % PROBFEM-based inference
    proposalPROB = @(w) GaussianProposal(w, sigmaPropProb);
%     likelihoodPROB = @(w) ProbTwoLikelihood(w, xObs, uObs, obsNoise, f, x, nMC);
    likelihoodPROB = @(w) ProbTwoLikelihood(w, xObs, uObs, obsNoise,  A, F, x, nMC);
    posteriorPROB = @(w) likelihoodPROB(w) + prior(w);
    thetaPROB = MetropolisHastings(priorMean, posteriorPROB, proposalPROB, nMCMC, true);
    thetaPROB = thetaPROB(floor(nMCMC / 10) : end);
    estDistProb = [estDistProb abs(mean(thetaPROB) - exEstimate)];
end

% loglog(1./NVec, estDist)
loglog(1./NVec, estDistProb)
hold on
loglog(1./NVec, 1./NVec.^2)

% % Plot distributions
% figure
% hold on
% histogram(theta, 'normalization', 'pdf')
% histogram(thetaFEM, 'normalization', 'pdf')
% title('FEM')
% figure
% hold on
% histogram(theta, 'normalization', 'pdf')
% histogram(thetaPROB, 'normalization', 'pdf')
% title('PROB')
