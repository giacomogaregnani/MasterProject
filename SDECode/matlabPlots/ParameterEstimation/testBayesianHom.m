% clc; clear; close all
%%

% x = dlmread('MultiMulti.txt');
% sol = dlmread('MultiMultiSol.txt');
x = dlmread('MultiHomo.txt');
sol = dlmread('MultiHomoSol.txt');

%%
T = 1;

nChains = 1;
hom = x(1, :);
trueVals = [1, 0.5];

sampleTot = x(2:end, 2:end);
sampleSize = size(sampleTot, 1) / nChains;
nParam = size(sampleTot, 2);

%%

if (size(sol, 1) < size(sol, 2))
    sol = sol';
end
t = linspace(0, T, size(sol, 1));
figure
plot(t, sol(:, 1:end-1));
hold on
plot(t, sol(:, end-1)+2*sol(:,end), 'k--');
plot(t, sol(:, end-1)-2*sol(:,end), 'k--');

% sampleTot = exp(sampleTot);
% trueVals = exp(trueVals);

aVals = linspace(min(sampleTot(:, 1)), max(sampleTot(:, 1)), 200);
sVals = linspace(min(sampleTot(:, 2)), max(sampleTot(:, 2)), 200);
[aa, ss] = meshgrid(aVals, sVals);
bwa = 2 * (4 * std(sampleTot(:, 1))^5 / (3 * size(sampleTot, 1)))^(1/5);
bws = 2 * (4 * std(sampleTot(:, 2))^5 / (3 * size(sampleTot, 1)))^(1/5);

F = mvksdensity(sampleTot, [aa(:), ss(:)], 'bandwidth', [bwa, bws]);

W = 6; H = 6;
fig = createFigure(W, H, 'enhanced', 2);
scatter(sampleTot(:, 1), sampleTot(:, 2), 10, '.')
hold on
contour(aa, ss, reshape(F, 200, 200));
hold on
scatter(log(trueVals(1)), log(trueVals(2)), 'xk', 'linewidth', 1)
xlim([-3, 2]);
ylim([-1.5, 0.5]);
xlabel('$\log(\alpha)$', 'interpreter', 'latex')
ylabel('$\log(\sigma)$', 'interpreter', 'latex')
% title({'Posterior $\mu^0$'}, 'interpreter', 'latex')
title({'Posterior $\mu^\varepsilon$'}, 'interpreter', 'latex')
% title({'Posterior $\tilde \mu^0$'}, 'interpreter', 'latex')
axis square
box on
% export_fig(fig, '../../../Reports/DraftMultiSDE_19/VERSION2/Figures/MultiHomo.eps', '-nocrop', '-painters')
% export_fig(fig, '../../../Reports/DraftMultiSDE_19/VERSION2/Figures/MultiMulti.eps', '-nocrop', '-painters')
% export_fig(fig, '../../../Reports/DraftMultiSDE_19/VERSION2/Figures/MultiHomoMod.eps', '-nocrop', '-painters')

% for j = 1 : nParam
%     figure
%     f(2*j-1) = axes;
%     hold on
%     
%     figure
%     f(2*j) = axes;
%     hold on
% end
% 
% for j = 1 : nParam
%     maxF = 0;
%     
%     sample = sampleTot(:, j);
%     sample = exp(sample);
%     
%     plot(f(nParam+j), sample);
%     
%     means = mean(sample);
%     stddevs = std(sample);
%     disp(['parameter ', num2str(j), ', mean = ' num2str(means), ', stddev = ', num2str(stddevs)])
%     
%     kspoints = linspace(means-3*stddevs, means+3*stddevs, 1000);
%     
%     if j == 3
%         histogram(f(j), sample)
%     else
%         [dens, points] = ksdensity(sample, kspoints);
%         plot(f(j), points, dens)
%     end
%     
%     idx = idx + sampleSize;
%     maxF = max(maxF, max(dens));
%     
%     plot(f(j), [hom(j), hom(j)], [0, maxF], 'k--');
%     plot(f(j), [trueVals(j), trueVals(j)], [0, maxF], 'k');
% end