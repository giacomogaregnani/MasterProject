clc; clear; close all;
%%
allResults = dlmread('convFEM.txt');

h = allResults(:, 1);
strongError = allResults(:, 2);
weakError = allResults(:, 3);

figure
loglog(h, strongError, 'o-')
hold on
loglog(h, weakError, 'o-')
loglog(h, h.^0.5 / 10, 'k')
loglog(h, h, 'k--')
legend('strong', 'weak', 'order 0.5', 'order 1', 'location', 'best')

strongOrder = log(strongError(1:end-1) ./ strongError(2:end)) / log(2);
weakOrder = log(weakError(1:end-1) ./ weakError(2:end)) / log(2);

title(['strong = ' num2str(mean(strongOrder)) ' ,weak = ' num2str(mean(weakOrder))]);