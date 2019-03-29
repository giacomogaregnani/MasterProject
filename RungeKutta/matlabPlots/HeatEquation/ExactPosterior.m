clc; clear; close all;

%% INVERSE PROBLEM PARABOLIC TEST

%% Data

dx = 0.01;
T = 0.01;
x = 0 : dx : 1;
N = length(x) - 1;

A = gallery('tridiag', N-1);
A = A / (dx^2);

u0 = (x < 0.9) .* (x > 0.7) + (x < 0.3) .* (x > 0.1);
% u0 = x .* (1 - x);
% u0 = ones(size(x));
u0 = u0';

%% Generate observations (gamma = 1)

% Forward operator of the heat equation
forwMatrix = expm(-A * T);

% Observation noise distribution
nEig = 1e4;
[V, D] = laplEig(x, nEig);
delta = 1e-8;
sqrtD = delta ./ diag(sqrt(-D));

noise = KL(zeros(size(x')), V, sqrtD, randn(nEig, 1));
figure
plot(x, noise)
title('noise realization')
yObs = [0; forwMatrix * u0(2:end-1); 0] + noise;

figure
plot(x, yObs)
title('observation')

%% Prior definition (beta = 1, alpha = 1, m0 = 0 => C_0 = -\Delta^(-1))

m0 = zeros(N+1, 1);

%% Posterior

SupportMatrix = eye(N-1) + 1 / delta * forwMatrix * forwMatrix;
m = m0(2:end-1) + 1 / delta * forwMatrix * (SupportMatrix \ (yObs(2:end-1) - forwMatrix * m0(2:end-1)));
m = [0; m; 0];
C = (A \ eye(N-1)) / SupportMatrix;

figure
hold on

% Draw samples from the posterior
S = chol(C);
for i = 1 : 1000
%     postSample = m + S' * randn(size(m));
    postSample = [0; m(2:end-1) + S' * randn(N-1, 1); 0];
    if i == 1
        p1 = plot(x, postSample, 'color', [0.6, 0.6, 0.6]);
    else
        plot(x, postSample, 'color', [0.6, 0.6, 0.6]);
    end
end

p2 = plot(x, m, 'k', 'linewidth', 2);
p3 = plot(x, u0, 'k--', 'linewidth', 2);

legend([p1, p2, p3], {'sample', 'mean', 'truth'})

