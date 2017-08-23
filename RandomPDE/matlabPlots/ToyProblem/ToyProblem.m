clc; clear; close all
%% Compare FEM and analytic solution for -u'' = theta

h = 0.2;
xGrid = 0 : h : 1;
nElements = length(xGrid) - 1;

A = 1 / h * gallery('tridiag', nElements - 1);
F = h * ones(nElements - 1, 1);

% Exact solution
thetaEx = 2;
uEx = @(x, theta) theta / 2 * (-x.^2 + x);

% Observation
gamma = 10^(-3);
xObs = 2 / 3;
uObs = uEx(xObs, thetaEx) + gamma * randn(1);
 

%% Inverse Problem

% Prior (N(mT, sT^2)) 
mThetaPrior = 10;
sThetaPrior = 1;
prior = @(theta) exp(-(theta - mThetaPrior)^2 / (2 * sThetaPrior^2));

% Modeling error
nModError = 1000;
thetaModError = random('normal', mThetaPrior, sThetaPrior, nModError, 1);
for i = 1 : nModError
   e(i) = uEx(xObs, thetaModError(i)) - interp1(xGrid, uEx(xGrid, thetaModError(i)), xObs);
end
modErrorMean = mean(e);
modErrorVar = var(e);

% Likelihood
likelihood = @(u) exp(-1 / (2 * gamma^2) * (uObs - u)^2);
likelihoodFEM = @(u) exp(-1 / (2 * (gamma^2 + modErrorVar)) * (uObs - u - modErrorMean)^2);

% MCMC
nMCMC = 2e4; nMC = 1;
known = true;
if known
    thetaAllFEM = Toy_MCMC(h, xGrid, prior, likelihood, xObs, nMCMC, nMC, 'FEM', uEx, true, h);
    thetaAllFEM = thetaAllFEM(floor(nMCMC/10):end); % Burn in
    thetaAllFEMCorr = Toy_MCMC(h, xGrid, prior, likelihoodFEM, xObs, nMCMC, nMC, 'FEM', uEx, true, h);
    thetaAllFEMCorr = thetaAllFEMCorr(floor(nMCMC/10):end); % Burn in
    thetaAllEx = Toy_MCMC(h, xGrid, prior, likelihood, xObs, nMCMC, nMC, 'EX', uEx, true, h);
    thetaAllEx = thetaAllEx(floor(nMCMC/10):end);
end
thetaAllProb = Toy_MCMC(h, xGrid, prior, likelihood, xObs, nMCMC, nMC, 'PROB', uEx, false, h);
thetaAllProb = thetaAllProb(floor(nMCMC/10):end);


%% Plot results
figure
hold on
if known
    histogram(thetaAllEx, 'normalization', 'pdf', 'edgecolor', 'none')
    histogram(thetaAllFEM, 'normalization', 'pdf', 'edgecolor', 'none')
    histogram(thetaAllFEMCorr, 'normalization', 'pdf', 'edgecolor', 'none')
end
histogram(thetaAllProb, 'normalization', 'pdf', 'edgecolor', 'none')

% Theoretical
% exact
uTildeEx = uEx(xObs, thetaEx) / thetaEx;
mu = (uTildeEx * uObs * sThetaPrior^2 + mThetaPrior * gamma^2) ...
    / (gamma^2 + uTildeEx^2 * sThetaPrior^2);
varEx = sThetaPrior^2 * gamma^2 / (gamma^2 + uTildeEx^2 * sThetaPrior^2);
xPlot = linspace(mu - 5 * sqrt(varEx), mu + 5 * sqrt(varEx), 1000);
exDistr = normpdf(xPlot, mu, sqrt(varEx));
plot(xPlot, exDistr, 'k--', 'linewidth', 2)

% FEM
uTilde = interp1(xGrid, uEx(xGrid, thetaEx), xObs) / thetaEx;
muFEM = (uTilde * uObs * sThetaPrior^2 + mThetaPrior * gamma^2) ...
    / (gamma^2 + uTilde^2 * sThetaPrior^2);
varFEM = sThetaPrior^2 * gamma^2 / (gamma^2 + uTilde^2 * sThetaPrior^2);
xPlotFEM = linspace(muFEM - 5 * sqrt(varFEM), muFEM + 5 * sqrt(varFEM), 1000);
plot(xPlotFEM, normpdf(xPlotFEM, muFEM, sqrt(varFEM)), 'k', 'linewidth', 2)

% FEM CORRECTED
muFEMCorr = ((gamma^2 + modErrorVar) * mThetaPrior + uTilde * (uObs - modErrorMean) * sThetaPrior^2) / ...
    (gamma^2 + modErrorVar + sThetaPrior^2 * uTilde^2);
varFEMCorr =  (sThetaPrior^2 * (gamma^2 + modErrorVar))/ ...
    (gamma^2 + modErrorVar + sThetaPrior^2 * uTilde^2);
xPlotFEMCorr = linspace(muFEMCorr - 5 * sqrt(varFEMCorr), muFEMCorr + 5 * sqrt(varFEMCorr), 1000);
plot(xPlotFEMCorr, normpdf(xPlotFEMCorr, muFEMCorr, sqrt(varFEMCorr)), 'k-.', 'linewidth', 2)

% PROB
% tMin = min(thetaAllProb); tMax = max(thetaAllProb);
% nTSample = 10000;
% tSample = linspace(tMin - 0.1, tMax + 0.1, nTSample);
% priorProb = exp(-(tSample - mThetaPrior).^2 / (2 * sThetaPrior^2));
% likProbIntegrand = @(h, theta) h.^(-2) .* exp(-0.5 * gamma^(-2) * (uObs - theta / 4 * h .* (1 - h)).^2);
% likProb = zeros(size(priorProb));
% hSample = linspace(h - h^(2) / 2, h + h^2 / 2, 10000);
% for i = 1 : nTSample
%    likProb(i) = trapz(hSample, likProbIntegrand(hSample, tSample(i)));
% end
% postProb = priorProb .* likProb;
% postProb = postProb / trapz(tSample, postProb);
% 
% plot(tSample, postProb, 'k:', 'LineWidth', 2);
% 
% set(gca, 'yTick', [])

% Truth
plot([thetaEx, thetaEx], [0, max(exDistr)], 'r--', 'LineWidth', 2)  

if known
    legend('EX-MCMC', 'FEM-MCMC', 'FEM-CORR', 'PROB-MCMC', 'EX-EX', 'FEM-EX', 'FEM-CORR-EX', 'TRUTH', 'location', 'best')
else
    legend('PROB-MCMC', 'EX-EX', 'FEM-EX', 'PROB-EX', 'TRUTH', 'location', 'best')
end