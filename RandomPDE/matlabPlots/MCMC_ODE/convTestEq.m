clc; clear; close all

%% Construct exact posterior

obs = dlmread('observations.txt');
tObs = obs(:, 1);
yObs = obs(:, 2);

prior = @(theta) exp(-theta^2 / 2);
sigmaObs = 0.005;
likelihood = @(theta) computeLikelihood(theta, tObs, yObs, sigmaObs);
posterior = @(theta) prior(theta) * likelihood(theta);


%% Prob FEM posteriors

figure
hold on

hSample = 0.001;

HellDet = zeros(9, 1);
for i = 1 : length(HellDet)
    x = dlmread(['MCMC_ODE_DET_4_TEST', num2str(i), '.txt']);
    
    thetaSample = min(-1.05, 1.1 * min(x)) : hSample : max(-0.95, 0.9 * max(x));
    postEx = zeros(size(thetaSample));
    nSample = length(thetaSample);
    
    for j = 1 : nSample
        postEx(j) = posterior(thetaSample(j));
    end
    postEx = postEx / trapz(thetaSample, postEx);
    
    if i == 1
        plot(thetaSample, postEx, 'k', 'linewidth', 2)
    end
        
    [f, ~] = ksdensity(x, thetaSample);
    plot(thetaSample, f)
       
    HellDet(i) = 1 / sqrt(2) * sqrt(hSample * sum((sqrt(f) - sqrt(postEx)).^2));
end
title('DET')

figure
hold on

HellProb = zeros(4, 1);
for i = 1 : length(HellProb)
    x = dlmread(['MCMC_ODE_PROB_5_TEST', num2str(i), '.txt']);
    
    thetaSample = min(-1.05, 1.1 * min(x)) : hSample : max(-0.95, 0.9 * max(x));
    nSample = length(thetaSample);
    postEx = zeros(size(thetaSample));
    
    for j = 1 : nSample
        postEx(j) = posterior(thetaSample(j));
    end
    postEx = postEx / trapz(thetaSample, postEx);
    
    if i == 1
        plot(thetaSample, postEx, 'k', 'linewidth', 2)
    end
    
    [f, ~] = ksdensity(x, thetaSample);
    plot(thetaSample, f)
    HellProb(i) = 1 / sqrt(2) * sqrt(hSample * sum((sqrt(f) - sqrt(postEx)).^2));
       
end
title('PROB')

%%

hDet = 0.5 * 2.^[0:-1:-length(HellDet)+1];
hProb = 0.5 * 2.^[0:-1:-length(HellProb)+1];
figure
loglog(hProb, HellProb, 'o-')
hold on
loglog(hDet, HellDet, 'o-')
loglog(hProb, 40 * hProb .^1.1, 'k')
loglog(hDet, 5 * hDet .^2, 'k--')
legend('PROB', 'DET', 'ORDER', 'location', 'best')