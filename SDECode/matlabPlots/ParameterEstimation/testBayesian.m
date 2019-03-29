clc; clear; close all
%%

x = dlmread('testHomo.txt');
sol = dlmread('testHomoSol.txt');

%%
T = 10;
t = linspace(0, T, size(sol, 2));
figure
plot(t, sol');
legend('true path', 'noisy path')

nChains = 4;
hom = x(1, :);
trueVals = [1, 0.5];

sampleTot = x(2:end, 2:end);
sampleSize = size(sampleTot, 1) / nChains;

figure
f(1) = axes;
hold on

% figure
% f(2) = axes;
% hold on

figure
f(3) = axes;
hold on

% figure
% f(4) = axes;
% hold on

for j = 1 : 1
    idx = 1;
    maxF = 0;
    
    for i = 1  : nChains
        
        sample = sampleTot(idx:idx+sampleSize-1, j);
%         if j == 21
            sample = exp(sample);
%         end
        plot(f(2+j), sample);        
        
        [dens, points] = ksdensity(sample);
%         histogram(f(j), sample, 'normalization', 'pdf')
        
        plot(f(j), points, dens)
        idx = idx + sampleSize;
        
        maxF = max(maxF, max(dens));
        
    end
    
    plot(f(j), [hom(j), hom(j)], [0, maxF], 'k--');
    plot(f(j), [trueVals(j), trueVals(j)], [0, maxF], 'k');
    
end