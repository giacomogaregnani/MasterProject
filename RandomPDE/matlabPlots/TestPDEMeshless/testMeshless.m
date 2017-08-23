clc; clear; close all;
%% Test of meshless method.

l = 0.6;

xObs = [0, rand(1, 5), 1];
xObsDomain = xObs(2:end-1);
xObsBound = [xObs(1), xObs(end)];

f = @(x) sin(2 * pi * x);
% f = @(x) x;
uEx = @(x) -1 / (2 * pi)^2 * sin(2 * pi * x);
% uEx = @(x) -1 / (2 * pi)^2 * (cos(2 * pi * x) - 1);

F = f(xObsDomain)';
B = [0; 0];

xEval = linspace(0, 1, 100);

LLbK = assembleLLbK(xObsDomain, xObsBound, l);
LK = assembleLK(xObsDomain, xObsBound, xEval, l);
LbK = assembleLbK(xObsDomain, xObsBound, xEval, l);
K = assembleK(xEval, xEval, l);

mu = LbK * (LLbK \ [F; B]);
Sigma = K - LbK * (LLbK \ LK);

plot(xEval, mu, 'b-');
hold on
plot(xEval, mu + 2 * sqrt(diag(Sigma)), 'b--')
plot(xEval, mu - 2 * sqrt(diag(Sigma)), 'b--')
plot(xEval, uEx(xEval), 'k')
plot(xObs, uEx(xObs), 'ro')