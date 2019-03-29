clc; clear; close all
%%
truth = dlmread('infDimtruth.txt');
resultProb = dlmread('infDimprob.txt');
resultDet = dlmread('infDim.txt');

%%
xTruth = linspace(0, 1, length(truth));

meanDet = resultDet(1, :);
stdDet = sqrt(diag(resultDet(2:end, :))');

meanProb = resultProb(1, :);
stdProb = sqrt(diag(resultProb(2:end,:))');

x = linspace(0, 1, length(meanDet));
x2 = [x, fliplr(x)];

figure
hold on
confIntProb = [meanProb + 1.96 * stdProb, fliplr(meanProb - 1.96 * stdProb)];
fill(x2, confIntProb, 0.9 * ones(1, 3), 'linestyle', 'none')
% plot(x, meanProb + 1.96 * stdProb, 'r--')
% plot(x, meanProb - 1.96 * stdProb, 'r--')

confIntDet = [meanDet - 1.96 * stdDet, fliplr(meanDet + 1.96 * stdDet)];
fill(x2, confIntDet, 0.7 * ones(1, 3),'linestyle', 'none')
plot(x, meanDet, 'k--')
plot(x, meanProb, 'k')
% hold on
% plot(x, meanDet + 1.96 * stdDet, 0.7 * ones(1, 3))
% plot(x, meanDet - 1.96 * stdDet, 0.7 * ones(1, 3))

plot(xTruth, truth, 'k', 'linewidth', 2)


% 
% x = 1 : 300;
% curve1 = log(x);
% curve2 = 2*log(x);
% plot(x, curve1, 'r', 'LineWidth', 2);
% hold on;
% plot(x, curve2, 'b', 'LineWidth', 2);
% x2 = [x, fliplr(x)];
% inBetween = [curve1, fliplr(curve2)];
% fill(x2, inBetween, 'g');

% ylim([0.6, 1.4])

%% Cut stuff
% 
% valuesTrue = exp([0, -0.3, 0.3]);
% values = [1, 0.8, 1.2];
% 
% meanDetCut = valuesTrue(2) * (meanDet < values(2)) ...
%            + valuesTrue(3) * (meanDet > values(3)) ...
%            + valuesTrue(1) * (meanDet <= values(3)) .* (meanDet >= values(2));
% 
% meanProbCut = valuesTrue(2) * (meanProb < values(2)) ...
%             + valuesTrue(3) * (meanProb > values(3)) ...
%             + valuesTrue(1) * (meanProb <= values(3)) .* (meanDet >= values(2));
%        
% figure
% plot(x, meanDetCut, 'b')
% hold on
% plot(x, meanProbCut, 'r')
% plot(xTruth, truth, 'k')

%% Read KL coeffs

% coeffsDet = dlmread('infDimcoeffs.txt');
% coeffsProb = dlmread('infDimcoeffsProb.txt');
% % coeffTrue = [1,1,0,1];
% 
% col = parula(size(coeffsDet,2)*2);
% col = col(1:2:end, :);
% for i = 1 : size(coeffsDet,2)
%    figure
%    hold on
%    [dens, xi] = ksdensity(coeffsProb(:, i));
%    plot(xi, dens, '--', 'color', col(i, :))
%    [dens, xi] = ksdensity(coeffsDet(:, i));
%    plot(xi, dens, 'color', col(i, :))
% %    plot([coeffTrue(i), coeffTrue(i)], [0, max(dens)], 'k--')
% end