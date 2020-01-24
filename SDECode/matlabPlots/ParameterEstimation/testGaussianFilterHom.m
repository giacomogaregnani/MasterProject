clc; clear; close all

%% CAREFUL ABOUT ALPHA
%%
results = dlmread('GaussianTest_SI_HOM.txt');
n = length(results);
x = results(1:n/2);
z = results(n/2+1:end);
clear results

% figure
% scatter(x, z, 10, '.')

% v = @(x) x.^2 / 2;
% v = @(x) x.^4 / 4;
v = @(x) x.^6 / 6;
% v = @(x) x.^4 / 4 - x.^2 / 2;
alpha = 3; sigma = 0.5;

figure
[f, xi] = ksdensity(x, 'NumPoints', 2000, 'bandwidth', 0.1);
plot(xi, f)
hold on
[f,xi] = ksdensity(z, 'NumPoints', 2000, 'bandwidth', 0.1);
plot(xi, f)

xwide = linspace(-5, 5, 10000);

vX = @(x) alpha / sigma * v(x);
theorypdf = @(x) exp(-vX(x));
inttheorypdf = trapz(xwide, theorypdf(xwide));
theorypdfnorm = @(x) theorypdf(x) / inttheorypdf;

vZ = @(z) alpha^3 * (1 + alpha) / sigma * v(z);
theorypdflimit = @(z) exp(-vZ(z));
inttheorypdflimit = trapz(xwide, theorypdflimit(xwide));
theorypdflimitnorm = @(z) theorypdflimit(z) / inttheorypdflimit;

plot(xi, theorypdfnorm(xi))
plot(xi, theorypdflimitnorm(xi))
% legend('result', 'theory')
legend('true', 'filter', 'theory', 'theory filter')

% theoryJoint = @(x, z) theorypdfnorm(x) .* theorypdflimitnorm(z) .* exp(-1/sigma*(x-z).^4);
% [XX, ZZ] = meshgrid(xi,xi);
% theoryJointEval = theoryJoint(XX, ZZ);
% 
% figure
% contour(XX, ZZ, theoryJointEval, 20)
% hold on
% scatter(x, z, 10, '.')