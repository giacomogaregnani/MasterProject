clc; clear; close all
%% Generate observations

thetaEx = 1;
f = @(x) sin(2 * pi * x)';
xRef = linspace(0, 1, 10000);
uEx = @(x) 1 / (2 * pi)^2 * sin(2 * pi * x);

% Generate observations
xObs = 0.8;
obsNoise = 1e-5;
uObs = 1 / exp(thetaEx) * uEx(xObs) + obsNoise * randn(1);

%% Inference

nMCMC = 2e4;
sigmaProp = 0.005;
proposal = @(w) GaussianProposal(w, sigmaProp);
priorMean = 0;
priorVariance = 1;
prior = @(w) GaussianPrior(w, priorMean, priorVariance);

% Exact inference
likelihoodEx = @(w) ExactLikelihood(w, xObs, uObs, obsNoise);
posteriorEx = @(w) likelihoodEx(w) + prior(w);
theta = MetropolisHastings(priorMean, posteriorEx, proposal, nMCMC);
theta = theta(floor(nMCMC / 10) : end);

% FEM-based inference
N = 4; h = 1 / N; x = linspace(0, 1, N+1);
A = 1 / h * gallery('tridiag', N-1);
F = h * f(x(2:end-1));
likelihoodFEM = @(w) FEMLikelihood(w, xObs, uObs, obsNoise, A, F, x);
posteriorFEM = @(w) prior(w) + likelihoodFEM(w);
thetaFEM = MetropolisHastings(priorMean, posteriorFEM, proposal, nMCMC);
thetaFEM = thetaFEM(floor(nMCMC / 10) : end);

% % PROBFEM-based inference
% nMC = 10;
% likelihoodPROB = @(w) ProbLikelihood(w, xObs, uObs, obsNoise, f, x, nMC);
% posteriorPROB = @(w) likelihoodPROB(w) + prior(w);
% thetaPROB = MetropolisHastings(priorMean, posteriorPROB, proposal, nMCMC);
% thetaPROB = thetaPROB(floor(nMCMC / 10) : end);

% % Kaipio-Somersalo method
% nTrials = 100;
% [m, sigmaErr] = estimateErrorStats(nTrials, xObs, A, F, x, priorMean, priorVariance, uEx);
% likelihoodKS = @(w) KSLikelihood(w, xObs, uObs, obsNoise, A, F, x, m, sigmaErr);
% posteriorKS = @(w) likelihoodKS(w) + prior(w);
% thetaKS = MetropolisHastings(priorMean, posteriorKS, proposal, nMCMC);
% % thetaKS = thetaKS(floor(nMCMC / 10) : end);
% 
% % Mixed method
% likelihoodMIX = @(w) ProbLikelihood(w, xObs, uObs - m, obsNoise, f, x, nMC);
% posteriorMIX = @(w) likelihoodMIX(w) + prior(w);
% thetaMIX = MetropolisHastings(priorMean, posteriorMIX, proposal, nMCMC);
% thetaMIX = thetaMIX(floor(nMCMC / 10) : end);

% Iterated Calvetti-Somersalo method
nIterIKS = 10;
thetaIKS = priorMean;
nTrials = 500;
for i = 1 : nIterIKS
    if i == 1
        [m, sigmaErr] = estimateErrorStats(nTrials, xObs, A, F, x, priorMean, priorVariance, uEx);
        display(['m = ', num2str(m), ' sigmaErr = ', num2str(sigmaErr)])
    else
%         prior = @(w) GaussianPrior(w, mean(thetaIKS), var(thetaIKS));
        [m, sigmaErr] = estimateErrorStats(nTrials, xObs, A, F, x, mean(thetaIKS), var(thetaIKS), uEx);
        display(['m = ', num2str(m), ' sigmaErr = ', num2str(sigmaErr)])
    end
    likelihoodKS = @(w) KSLikelihood(w, xObs, uObs, obsNoise, A, F, x, m, sigmaErr);
    posteriorKS = @(w) likelihoodKS(w) + prior(w);
    thetaIKS = [thetaIKS, MetropolisHastings(thetaIKS(end), posteriorKS, proposal, nMCMC / nIterIKS)];
end
thetaIKS = [thetaIKS, MetropolisHastings(thetaIKS(end), posteriorKS, proposal, nMCMC)];
thetaIKS = thetaIKS(floor(length(thetaIKS) / 10) : end);

% Prob Variant Iterated Calvetti-Somersalo method
display('prob')
thetaIKSProb = priorMean;
nMC = 10;
nTrials = 500;
for i = 1 : nIterIKS
    if i == 1
        [mProb, sigmaErrProb] = estimateErrorStatsProb(nTrials, xObs, x, f, priorMean, priorVariance, uEx, nMC);
        display(['m = ', num2str(mProb), ' sigmaErr = ', num2str(sigmaErrProb)])
    else
%         prior = @(w) GaussianPrior(w, mean(thetaIKSProb), var(thetaIKSProb));
        [mProb, sigmaErrProb] = estimateErrorStatsProb(nTrials, xObs, x, f, mean(thetaIKSProb), var(thetaIKSProb), uEx, nMC);
        display(['m = ', num2str(mProb), ' sigmaErr = ', num2str(sigmaErrProb)])
    end
    likelihoodKS = @(w) KSLikelihoodProb(w, xObs, uObs, obsNoise, A, F, x, mProb, sigmaErrProb);
    posteriorKS = @(w) likelihoodKS(w) + prior(w);
    thetaIKSProb = [thetaIKSProb, MetropolisHastings(thetaIKSProb(end), posteriorKS, proposal, nMCMC / nIterIKS)];
end
thetaIKSProb = [thetaIKSProb, MetropolisHastings(thetaIKSProb(end), posteriorKS, proposal, nMCMC)];
thetaIKSProb = thetaIKSProb(floor(length(thetaIKSProb) / 10) : end);

% Plot distributions
figure
hold on
histogram(theta, 'normalization', 'pdf')
histogram(thetaFEM, 'normalization', 'pdf')
title('FEM')
% figure
% hold on
% histogram(theta, 'normalization', 'pdf')
% histogram(thetaPROB, 'normalization', 'pdf')
% title('PROB')
% figure
% hold on
% histogram(theta, 'normalization', 'pdf')
% histogram(thetaKS, 'normalization', 'pdf')
% title('KS')
% figure
% hold on
% histogram(theta, 'normalization', 'pdf')
% histogram(thetaMIX, 'normalization', 'pdf')
% title('MIX')
figure
hold on
histogram(theta, 'normalization', 'pdf')
histogram(thetaIKS, 'normalization', 'pdf')
title('IKS')
figure
hold on
histogram(theta, 'normalization', 'pdf')
histogram(thetaIKSProb, 'normalization', 'pdf')
title('IKSProb')

% Trace plots
% figure
% plot(theta)
% hold on
% plot(thetaFEM)
% plot(thetaKS)
% plot(thetaPROB)
% plot(thetaMIX)
% plot(thetaEx * ones(size(theta)), 'k--', 'linewidth', 2)