% Load
clear; clc; close all;


h = 0.05 * 2 .^ [0 : -1 : -11];

for i = 1 : length(h)
    
    EE = dlmread(['invMeasEE_h', num2str(floor(h(i)*1e6)), '.txt']);
    EI = dlmread(['invMeasEI_h', num2str(floor(h(i)*1e6)), '.txt']);
    MP = dlmread(['invMeasMP_h', num2str(floor(h(i)*1e6)), '.txt']);
    GA = dlmread(['invMeasGA_h', num2str(floor(h(i)*1e6)), '.txt']);
    nComponents = size(EE, 2);
    
    for j = 1 : nComponents
        figure
        subplot(2, 2, 1)
        histogram(EE(:, j), 20, 'Normalization', 'pdf')
        subplot(2, 2, 2)
        histogram(EI(:, j), 20, 'Normalization', 'pdf')
        subplot(2, 2, 3)
        histogram(MP(:, j), 20, 'Normalization', 'pdf')
        subplot(2, 2, 4)
        histogram(GA(:, j), 20, 'Normalization', 'pdf')
        
        varEE(i, j) = var(EE(:, j));
        varEI(i, j) = var(EI(:, j));
        varMP(i, j) = var(MP(:, j));
        varGA(i, j) = var(GA(:, j));
    end
    
end

for j = 1 : nComponents
    figure
    loglog(h, varEE(:, j), 'o-')
    hold on
    loglog(h, varEI(:, j),'o-')
    loglog(h, varMP(:, j),'o-')
    loglog(h, varGA(:, j),'o-')
    loglog(h, h.^2, 'k--')
    loglog(h, h.^4, 'k')
    loglog(h, h.^6, 'k--')
    loglog(h, h.^8, 'k')
    legend('EE', 'EI', 'MP', 'GA', '2', '4', '8', 'location', 'SE')
end