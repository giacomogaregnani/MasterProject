% clc; clear; close all

%%
results = dlmread('GaussianTest2.txt');
n = length(results);

x = results(1:n/2);
z = results(n/2+1:end);

x1 = -3:.02:3;
x2 = -3:.02:3;
[X1,X2] = meshgrid(x1,x2);

% F = mvksdensity([x, z], [X1(:), X2(:)], 'bandwidth', std([x, z])*(4/3/numel([x, z]))^(1/5));
% F = reshape(F, length(x1), length(x1));

figure
scatter(x, z, 10, '.')

% mu = [mean(x) mean(z)];
% sigma = cov(x, z);

alpha = 1; sigma = 0.5; delta = 0.05;
muTh = [0, 0];
C = 1 / (1 + alpha*delta);
SigmaTh = sigma/alpha * [1, C; C, C];

X = [X1(:) X2(:)];
Y = mvnpdf(X, muTh, SigmaTh);
Y = reshape(Y, length(x2), length(x1));

hold on
% contour(x1, x2, Y, 20)
xlabel('x')
ylabel('z')
% contour(x1, x2, F, 10)

figure
[f,xi] = ksdensity(x, 'bandwidth', 0.2);
plot(xi, f)
hold on
[f,xi] = ksdensity(z, 'bandwidth', 0.2);
plot(xi, f)
% 
% p = @(x) sin(x);
% v = @(x) x.^2 / 2;
% alpha = 1;
% sigma = 0.5;
% eps = 0.1;
% 
% xwide = linspace(-10, 10, 10000);
% 
% theorypdf = @(x) exp(-alpha / sigma * v(x) - 1 / sigma * p(x/eps));
% theorypdfnorm = @(x) theorypdf(x) / trapz(xwide, theorypdf(xwide));
% 
% plot(xi, theorypdfnorm(xi))
% 
% legend('true', 'filter', 'theory')