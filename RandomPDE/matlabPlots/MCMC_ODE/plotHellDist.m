clc; clear; close all
%%

filenameRef = 'MCMC_ODE_HELL_REF_h2';
filename = 'MCMC_ODE_HELL_NOISY_h2';

xRef = dlmread([filenameRef, '.txt']);

nSamplesKD = 1e3;
sample = linspace(min(xRef), max(xRef), nSamplesKD);
h = (sample(end) - sample(1)) / nSamplesKD;
[fRef, xi] = ksdensity(xRef, sample);

nMC_H = 10;
H = zeros(1, 7);

for i = 1 : length(H)
    display(['======== i = ', num2str(i)])

    hMC_H = zeros(1, nMC_H);
    
    parfor j = 0 : nMC_H-1
        
        display(['j = ' num2str(j)])
        x = dlmread([filename, num2str(i), '_', num2str(j), '.txt']);
        [f, xi] = ksdensity(x, sample);
        
        HCurr = 1 / sqrt(2) * sqrt(h * sum((sqrt(f) - sqrt(fRef)).^2));
        hMC_H(j+1) = HCurr;
        
    end
    
    H(i) = mean(hMC_H);

end



%%
nMC = 2.^[0:6];
figure
loglog(nMC, H, 'o-')
hold on
loglog(nMC, 1./sqrt(nMC))