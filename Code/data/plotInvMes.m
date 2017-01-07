% Load
clear; clc; close all;

h = 0.025 * 2 .^ [0 : -1 : -3];

for i = 1 : length(h)
    
    EE = dlmread(['invMeasEE_h', num2str(floor(h(i)*1e6)), '.txt']);
    EI = dlmread(['invMeasEI_h', num2str(floor(h(i)*1e6)), '.txt']);
    MP = dlmread(['invMeasMP_h', num2str(floor(h(i)*1e6)), '.txt']);
    GA = dlmread(['invMeasGA_h', num2str(floor(h(i)*1e6)), '.txt']);
    nComponents = size(EE, 2);
    
%     for j = 1 : nComponents
%         figure
%         subplot(2, 2, 1)
%         histogram(EE(:, j), 20, 'Normalization', 'pdf')
%         title('EE')
%         subplot(2, 2, 2)
%         histogram(EI(:, j), 20, 'Normalization', 'pdf')
%         title('EI')
%         subplot(2, 2, 3)
%         histogram(MP(:, j), 20, 'Normalization', 'pdf')
%         title('MP')
%         subplot(2, 2, 4)
%         histogram(GA(:, j), 20, 'Normalization', 'pdf')
%         title('GA')
%         
%         varEE(i, j) = var(EE(:, j));
%         varEI(i, j) = var(EI(:, j));
%         varMP(i, j) = var(MP(:, j));
%         varGA(i, j) = var(GA(:, j));
%     end
    
    figure
    hold on
    [f, xi] = ksdensity(EE(:, 1));
    plot(xi, f);
    [f, xi] = ksdensity(MP(:, 1));
    plot(xi, f);
    legend('EE', 'MP')
    figure
    hold on
    [f, xi] = ksdensity(EI(:, 1));
    plot(xi, f);
    [f, xi] = ksdensity(GA(:, 1));
    plot(xi, f);
    legend('EI', 'GA')
    
end