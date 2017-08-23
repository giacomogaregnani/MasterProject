clc; close all; clear
%% Generate data

f = @(x) exp(-sin(3*x).^2 - x.^2);
sigma = 0.05;

X = linspace(-3, 3, 15);
F = f(X)';
F = F + sigma * randn(size(X'));


%% Prior distribution

pF = struct('mu', [], 'k', []);
pF.mu = 0;
% pF.k = @(a, b) min(a, b);
pF.k = @(a, b) exp(-1 / (2 * 0.2) * abs(a - b).^2); 

%% Test points

Xstar = linspace(-5, 5, 500);

%% Plot prior

mu = zeros(size(Xstar))';
Sigma = buildKmatrix(Xstar, Xstar, pF);

plot(Xstar, mu, 'b', 'linewidth', 2);
hold on
plot(Xstar, mu + 2 * sqrt(diag(Sigma)), 'b--')
plot(Xstar, mu - 2 * sqrt(diag(Sigma)), 'b--') 


%% Matrices

Kxx = buildKmatrix(X, X, pF);
KxxStar = buildKmatrix(Xstar, X, pF);
KxStarxStar = buildKmatrix(Xstar, Xstar, pF);

%% Posterior parameters

mu = KxxStar * ((Kxx + sigma^2 * eye(size(Kxx))) \ F);
Sigma = KxStarxStar - KxxStar * ((Kxx + sigma^2 * eye(size(Kxx))) \ KxxStar');

 %% Plot posterior

% fig = createFigure(16, 7);
figure
plot(X, F, 'or', 'markersize', 8)
hold on
plot(Xstar, mu, 'linewidth', 2)
plot(Xstar, mu + 2 * sqrt(diag(Sigma)), 'b--')
plot(Xstar, mu - 2 * sqrt(diag(Sigma)), 'b--')
% true function
plot(Xstar, f(Xstar), 'k', 'linewidth', 2)
% 
% legend('observations', '\mu', '\mu + 2\Sigma_{ii}^{0.5}', ...
%        '\mu - 2\Sigma_{ii}^{0.5}', 'f', 'location', 'best')
% print -depsc2 ../../../Reports/Dobbiaco2017/posterior.eps

