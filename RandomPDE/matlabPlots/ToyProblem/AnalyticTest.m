clc; clear; close all
%% 
h = 0.1;
xGrid = 0 : h : 1;
nElements = length(xGrid) - 1;

A = 1 / h * gallery('tridiag', nElements - 1);
F = h * ones(nElements - 1, 1);

% Exact solution
thetaEx = 2;
uEx = @(x, theta) theta / 2 * (-x.^2 + x);

% Observation
gamma = 10^(-4);
xObs = 2 / 3;
uObs = uEx(xObs, thetaEx) + gamma * randn(1);
 
% Prior (N(mT, sT^2)) 
mThetaPrior = 2;
sThetaPrior = 1;
prior = @(theta) exp(-(theta - mThetaPrior)^2 / (2 * sThetaPrior^2));

% Likelihood
likelihood = @(u) exp(-1 / (2 * gamma^2) * (u - uObs)^2);

%% Plot results
figure
hold on

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

% Estimate modeling error mean and variance
nModError = 1000;
thetaModError = random('normal', mThetaPrior, sThetaPrior, nModError, 1);
for i = 1 : nModError
   e(i) = uEx(xObs, thetaModError(i)) - interp1(xGrid, uEx(xGrid, thetaModError(i)), xObs);
end
modErrorMean = mean(e);
modErrorVar = var(e);

muFEMCorr = ((gamma^2 + modErrorVar) * mThetaPrior + uTilde * (uObs - modErrorMean) * sThetaPrior^2) / ...
    (gamma^2 + modErrorVar + sThetaPrior^2 * uTilde^2);
varFEMCorr =  (sThetaPrior^2 * (gamma^2 + modErrorVar))/ ...
    (gamma^2 + modErrorVar + sThetaPrior^2 * uTilde^2);
xPlotFEMCorr = linspace(muFEMCorr - 5 * sqrt(varFEMCorr), muFEMCorr + 5 * sqrt(varFEMCorr), 1000);
plot(xPlotFEMCorr, normpdf(xPlotFEMCorr, muFEMCorr, sqrt(varFEMCorr)), 'k-.', 'linewidth', 2)

set(gca, 'yTick', [])

% Truth
plot([thetaEx, thetaEx], [0, max(exDistr)], 'r--', 'LineWidth', 2)  

legend('Ex', 'FEM', 'CORR', 'TRUTH')