% clc; clear; close all
%%

x = dlmread('MultiMulti.txt');
sol = dlmread('MultiMultiSol.txt');
% x = dlmread('MultiHomo.txt');
% sol = dlmread('MultiHomoSol.txt');
% fil = dlmread('testHomoFil.txt');

%%
T = 0.1;

nChains = 1;
hom = x(1, :);
trueVals = [1, 0.5];

sampleTot = x(2:end, 2:end);
sampleSize = size(sampleTot, 1) / nChains;
nParam = size(sampleTot, 2);

%%

t = linspace(0, T, size(sol, 2));
figure
plot(t, sol');
% figure
% plot(t, sol(end-1, :), 'k')
% hold on
% plot(t, sol(end-1, :) + 2*sol(end, :), 'k--')
% plot(t, sol(end-1, :) - 2*sol(end, :), 'k--')
% legend('error', 'IC')

for j = 1 : nParam
    figure
    f(2*j-1) = axes;
    hold on
    
    figure
    f(2*j) = axes;
    hold on
end

for j = 1 : nParam
    idx = 1;
    maxF = 0;
    
    for i = 1  : nChains
        
        sample = sampleTot(idx:idx+sampleSize-1, j);
        sample = exp(sample);
        
        
        plot(f(nParam+j), sample);
        
        means = mean(sample);
        stddevs = std(sample);
        disp(['parameter ', num2str(j), ', mean = ' num2str(means), ', stddev = ', num2str(stddevs)])
        
        kspoints = linspace(means-3*stddevs, means+3*stddevs, 1000);
        
        if j == 3
            histogram(f(j), sample)
        else
            [dens, points] = ksdensity(sample, kspoints);
            plot(f(j), points, dens)
        end
        
        idx = idx + sampleSize;       
        maxF = max(maxF, max(dens));
        
    end
    
    plot(f(j), [hom(j), hom(j)], [0, maxF], 'k--');
    plot(f(j), [trueVals(j), trueVals(j)], [0, maxF], 'k');
    
end