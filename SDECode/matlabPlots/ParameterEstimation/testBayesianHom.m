clc; clear; close all
%%

% x = dlmread('MultiMulti.txt');
% sol = dlmread('MultiMultiSol.txt');
x = dlmread('MultiHomo.txt');
sol = dlmread('MultiHomoSol.txt');
% x = dlmread('MultiHomoMod.txt');
% sol = dlmread('MultiHomoModSol.txt');

%%
T = 50;

nChains = 1;
hom = x(1, :);
trueVals = [1, 0.5];
% trueVals = [0.19, 0.09];

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
% hold on
% plot(t, sol(:, end-1)+2*sol(:,end), 'k--');
% plot(t, sol(:, end-1)-2*sol(:,end), 'k--');

% sampleTot = exp(sampleTot);
% trueVals = exp(trueVals);

aVals = linspace(min(sampleTot(:, 1))-0.5, max(sampleTot(:, 1))+0.5, 100);
sVals = linspace(min(sampleTot(:, 2))-0.2, max(sampleTot(:, 2))+0.2, 100);
[aa, ss] = meshgrid(aVals, sVals);

idx = 1:size(sampleTot, 1);
bwa = 2 * (4 * std(sampleTot(idx, 1))^5 / (3 * size(sampleTot, 1)))^(1/5);
bws = 2 * (4 * std(sampleTot(idx, 2))^5 / (3 * size(sampleTot, 1)))^(1/5);

F = mvksdensity(sampleTot(idx, :), [aa(:), ss(:)], 'bandwidth', [bwa, bws]);

% W = 6; H = 6;
% fig = createFigure(W, H, 'enhanced', 2);
figure
% scatter(sampleTot(idx, 1), sampleTot(idx, 2), 10, '.')
% hold on
contour(aa, ss, reshape(F, sqrt(length(F)), sqrt(length(F))));
hold on
scatter(log(trueVals(1)), log(trueVals(2)), 'xk', 'linewidth', 3)
% xlim([-3, 2]);
% ylim([-1.5, 0.5]);
xlabel('$\log(\alpha)$', 'interpreter', 'latex')
ylabel('$\log(\sigma)$', 'interpreter', 'latex')
% title({'Posterior $\mu^0$'}, 'interpreter', 'latex')
title({'Posterior $\mu^\varepsilon$'}, 'interpreter', 'latex')
% title({'Posterior $\tilde \mu^0$'}, 'interpreter', 'latex')
% axis square
% axis equal
box on
% export_fig(fig, '../../../Reports/DraftMultiSDE_19/VERSION2/Figures/MultiHomo.eps', '-nocrop', '-painters')
% export_fig(fig, '../../../Reports/DraftMultiSDE_19/VERSION2/Figures/MultiMulti.eps', '-nocrop', '-painters')
% export_fig(fig, '../../../Reports/DraftMultiSDE_19/VERSION2/Figures/MultiHomoMod.eps', '-nocrop', '-painters')

figure
plot(sampleTot(:, 1))
hold on
plot(log(trueVals(1)) .* ones(size(sampleTot(:, 1))), 'k')

figure
plot(sampleTot(:, 2))
hold on
plot(log(trueVals(2)) .* ones(size(sampleTot(:, 2))), 'k')

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
%     sample = sampleTot(idx, j);
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
%     maxF = max(maxF, max(dens));
%     
%     plot(f(j), [hom(j), hom(j)], [0, maxF], 'k--');
%     plot(f(j), [trueVals(j), trueVals(j)], [0, maxF], 'k');
% end