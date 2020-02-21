clc; clear; close all

%% CAREFUL ABOUT ALPHA
%%
results = dlmread('GaussianTest_SI_HOM.txt');
n = length(results);
x = results(1:n/2);
z = results(n/2+1:end);
clear results

figure
scatter(x, z, 10, '.')

% p = @(x) sin(x);
% v = @(x) x.^2 / 2;
p = @(x) cos(x);
v = @(x) x.^6 / 6;
% v = @(x) x.^4 / 4;
% v = @(x) x.^4 / 4 - x.^2 / 2;
alpha = 3; sigma = 0.5; eps = 0.05;

figure
[f, xi] = ksdensity(x, 'NumPoints', 2000, 'bandwidth', 0.1);
plot(xi, f)
hold on
[f,xi] = ksdensity(z, 'NumPoints', 2000, 'bandwidth', 0.1);
plot(xi, f)

xwide = linspace(-5, 5, 10000);

vEps = @(x) alpha / sigma * v(x) + 1 / sigma * p(x/eps);
theorypdf = @(x) exp(-vEps(x));
inttheorypdf = trapz(xwide, theorypdf(xwide));
theorypdfnorm = @(x) theorypdf(x) / inttheorypdf;

vEpsZ = @(z) alpha * (1 + 0.57) / sigma * v(z) + eps^2 / sigma * p(z / eps);
theorypdflimit = @(z) exp(-vEpsZ(z));
inttheorypdflimit = trapz(xwide, theorypdflimit(xwide));
theorypdflimitnorm = @(z) theorypdflimit(z) / inttheorypdflimit;

% plot(xi, theorypdfnorm(xi))
% plot(xi, theorypdflimitnorm(xi))
% legend('result', 'theory')
legend('true', 'filter', 'theory', 'theory filter')

% theoryJoint = @(x, z) theorypdfnorm(x) .* theorypdflimitnorm(z).* exp(-1/sigma*(alpha*x-z).^2);
% [XX, ZZ] = meshgrid(xi,xi);
% theoryJointEval = theoryJoint(XX, ZZ);
% 
% figure
% contour(XX, ZZ, theoryJointEval, 20)
% hold on
% scatter(x, z, 10, '.')