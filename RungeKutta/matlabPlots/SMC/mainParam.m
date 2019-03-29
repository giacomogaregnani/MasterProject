clc; clear; close all;
%% Define the ODE model and generate observations

exSol = [];

% Fitzhug Nagumo
% thetaEx = log([0.2, 0.2, 3.0]);
% y0 = [-1; 1];
% fParam = @(t, y, theta) [exp(theta(3)) * (y(1) - y(1)^3/3 + y(2));
%                          -1 / exp(theta(3)) * (y(1) - exp(theta(1)) + exp(theta(2)) * y(2))];

% Test equation
thetaEx = 1;
y0 = 1;
fParam = @(t, y, theta) theta(1) * y;
exSol = @(t, theta) exp(theta(1) * t);

fExParam = @(t, y) fParam(t, y, thetaEx);
tObs = 0.01:0.01:1;

if isempty(exSol)
    yRef = y0;
    tObs = [0, tObs];
    hRef = 1e-4;
    for k = 1 : length(tObs) - 1
        nSteps = 1 / hRef * (tObs(k+1) - tObs(k));
        yRef = [yRef rk4(fExParam, nSteps, yRef(:, end), hRef)];
    end
    tObs = tObs(2:end);
    yRef = yRef(:, 2:end);
else
    yRef = exSol(tObs, thetaEx);
end

% Corrupt data with noise
sigmaObs = 1;
yObs = yRef + sigmaObs * randn(size(yRef));

%% SMC
RK = @(f, N, y, h) explicitEuler(f, N, y, h);
RK2 = @(f, N, y, h) explicitEulerProb(f, N, y, h);
% RK = @(f, N, y, h) rkc(f, N, y, h, 10);
% RK2 = @(f, N, y, h) rkcProb(f, N, y, h, 10);
% RK2 = @(f, N, y, h) rk4(f, N, y, h);
obsProb = @(x, m) gaussianProb(x, sigmaObs, m);

nPar = 5e4;
h = 0.001;

[sample, theta, weights] = SMC_ODE_PAR_PROB(nPar, y0, fParam, RK, RK2, 0, tObs, yObs, h, obsProb, 1, length(thetaEx));
% [sample, theta, weights] = SMC_ODE_PAR(nPar, y0, fParam, RK, RK2, 0, tObs, yObs, h, obsProb, 0.97, length(thetaEx));

%% Plot solution histogram

relevant = find(weights > 0);

for i = 1 : length(y0)
    figure
    [f, xi] = ksdensity(sample(relevant, i), 'weights', weights(relevant));
    plot(xi, f)
    hold on
    plot([exSol(tObs(end), thetaEx(1)) exSol(tObs(end), thetaEx(1))], [0, max(f)], 'r--')
end
xlabel('state')

for i = 1 : length(thetaEx)
    figure
    [f, xi] = ksdensity(theta(relevant, i), 'weights', weights(relevant));
    plot(xi, f)
    hold on
    plot([thetaEx(i), thetaEx(i)], [0, max(f)], 'r--')
end
xlabel('parameter')