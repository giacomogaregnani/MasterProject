clc; clear; close all

%% Construct exact posterior

obs = dlmread('observations.txt');
tObs = obs(:, 1);
yObs = obs(:, 2);

prior = @(theta) exp(-theta^2 / 2);
sigmaObs = 0.005;
likelihood = @(theta) computeLikelihood(theta, tObs, yObs, sigmaObs);
posterior = @(theta) prior(theta) * likelihood(theta);

%%

nExp = 5;
h = 0.25 .* 2 .^ [0 : -1 : -(nExp-1)];
HellDet = zeros(nExp, 1);
HellProb = HellDet;

thetaSample = -3 : 0.001 : 0;
nSample = length(thetaSample);

postEx = zeros(1, nSample);
postET = postEx;
postETProb = postET;
p = 1;

figure
axDet = axes;
hold on
figure
axProb = axes;
hold on

for j = 1 : nSample
    postEx(j) = posterior(thetaSample(j));
end
postEx = postEx / trapz(thetaSample, postEx);
plot(axDet, thetaSample, postEx, 'k', 'linewidth', 2)
plot(axProb, thetaSample, postEx, 'k', 'linewidth', 2)

nMC = 1e3;
for i = 1 : nExp
    
    display(num2str(i))
    
    likelihoodEE = @(theta) computeLikelihoodEE(theta, tObs, yObs, sigmaObs, h(i));
    posteriorEE = @(theta) prior(theta) * likelihoodEE(theta);
    
    likelihoodEEProb = @(theta) computeLikelihoodEEProb(theta, tObs, yObs, sigmaObs, h(i), p, nMC);
    posteriorEEProb = @(theta) prior(theta) * likelihoodEEProb(theta);
    
    for j = 1 : nSample
        postET(j) = posteriorEE(thetaSample(j));
        postETProb(j) = posteriorEEProb(thetaSample(j));
    end
    postET = postET / trapz(thetaSample, postET);
    postETProb = postETProb / trapz(thetaSample, postETProb);
    
    plot(axDet, thetaSample, postET)
    plot(axProb, thetaSample, postETProb)
    
    HellDet(i) = 1 / sqrt(2) * sqrt(0.0001 * sum((sqrt(postET) - sqrt(postEx)).^2));
    HellProb(i) = 1 / sqrt(2) * sqrt(0.0001 * sum((sqrt(postETProb) - sqrt(postEx)).^2));
end

%%
figure
loglog(h, HellDet, 'o-')
hold on
loglog(h, HellProb, 'o-')
loglog(h, h, 'k--')
loglog(h, 50 * (h.^(2*p-1) + h.^(p-0.5)/nMC), 'k')
% loglog(h, h.^(2*p-1), 'k')