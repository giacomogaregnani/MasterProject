clc; clear; close all

load('homogeneous_T1.mat')

figure
plot(t, sol');
legend('true path', 'noisy path')

figure
hold on
idx = 1;
maxF = 0;

for i = 1  : nChains
    
    sample = sampleTot(idx:idx+sampleSize-1, 1);
    sample = exp(sample);
    
    [dens, points] = ksdensity(sample);
    
    plot(points, dens)
    idx = idx + sampleSize;
    
    maxF = max(maxF, max(dens));
    
end

plot([hom(1), hom(1)], [0, maxF], 'k--');
plot([trueVals(1), trueVals(1)], [0, maxF], 'k');
xlim = get(gca, 'xLim');

load('multiscale_T1.mat')

figure
plot(t, sol');
legend('true path', 'noisy path')

figure
hold on
idx = 1;
maxF = 0;

for i = 1  : nChains
    
    sample = sampleTot(idx:idx+sampleSize-1, 1);
    sample = exp(sample);
    
    [dens, points] = ksdensity(sample);
    
    plot(points, dens)
    idx = idx + sampleSize;
    
    maxF = max(maxF, max(dens));
    
end

plot([hom(1), hom(1)], [0, maxF], 'k--');
plot([trueVals(1), trueVals(1)], [0, maxF], 'k');
set(gca, 'xLim', xlim);

