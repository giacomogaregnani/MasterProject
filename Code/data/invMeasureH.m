clc; clear; close all;


% Variable h, deterministic
yLorenz = dlmread('invMeasEE_hTest.txt');
phi = @(y) sum(y.^2, 2);
phiY = phi(yLorenz);
h = 0.025 * (0.999 .^ [0 : 4999]);
intPhi = trapz(fliplr(h), fliplr(phiY));
intPhiNorm = 1 / (max(h)-min(h)) * intPhi;

phiYSorted = sort(phiY);
rangePhi = [0.5 * min(phiYSorted), 1.1 * max(phiYSorted)];
evalPoints = linspace(rangePhi(1), rangePhi(2), 1000); 
[f, xi] = ksdensity(phiY, evalPoints);
plot(xi, f, 'LineWidth', 2)
hold on
histogram(phiY, 'normalization', 'pdf', 'facealpha' , .5)
xlabel('\Phi(Y)')
ylabel('pdf')
title('deterministic, variable $h$')

intPhiMeasure = trapz(evalPoints, f .* evalPoints);

% Fixed h, probabilistic
yLorenzSto = dlmread('invMeasEE_hTestSto.txt');
phiYSto = phi(yLorenzSto);
intPhiStoNorm = mean(phiYSto);

phiYSorted = sort(phiYSto);
rangePhi = [0.5 * min(phiYSorted), 1.1 * max(phiYSorted)];
evalPoints = linspace(rangePhi(1), rangePhi(2), 10000); 
[f, xi] = ksdensity(phiYSto, evalPoints);
figure
plot(xi, f, 'LineWidth', 2)
hold on
histogram(phiYSto, 'normalization', 'pdf', 'facealpha' , .5)
xlabel('\Phi(Y)')
ylabel('pdf')
title('probabilistic, $M$ realizations')

intPhiMeasureSto = trapz(evalPoints, f .* evalPoints);