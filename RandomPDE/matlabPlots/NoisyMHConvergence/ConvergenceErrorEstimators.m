clc; clear; close all;
%% 

thetaEx = 1;
f = @(x) sin(2 * pi * x)';
xRef = linspace(0, 1, 10000);
uEx = @(x) 1 / (2 * pi)^2 * sin(2 * pi * x);
N = 5; h = 1 / N; x = linspace(0, 1, N+1);
A = 1 / h * gallery('tridiag', N-1);
F = h * f(x(2:end-1));

% Generate observations
xObs = 0.3;
obsNoise = 1e-5;
uObs = 1 / exp(thetaEx) * uEx(xObs) + obsNoise * randn(1);

% Specify prior
priorMean = 0;
priorVariance = 1;

% Compute "deterministic" error estimators
nTrials = 10000;
[m, sigmaErr] = estimateErrorStats(nTrials, xObs, A, F, x, priorMean, priorVariance, uEx);

nMC = [4, 8, 16, 32, 64, 128, 256, 512];
diff = [];

for n = nMC
    display(num2str(n));
    [mProb, sigmaErrProb] = estimateErrorStatsProb(n, xObs, x, f, priorMean, priorVariance, uEx, 10);
    diff = [diff, abs(m - mProb)];
end

loglog(nMC, diff)
