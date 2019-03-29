clc; clear; close all

filenameobs = {'observations1.txt', 'observations05.txt', ...
    'observations025.txt', 'observations0125.txt', 'observations0001.txt'};
enhanced = 1;
fontsizeLAB = getLatexTextSize('normalsize', 'enhanced', enhanced);
fontsizeTICK = getLatexTextSize('small', 'enhanced', enhanced);
W = 12; H = 9;

sigmaErr = [0.1, 0.05, 0.025, 0.0125, 0.0001];

col = colormap(parula(20));
close
col = col(1:4:end, :);

%%
fig = createFigure(W, H, 'enhanced',enhanced);
box on
hold on

hp = (0.5)^(1.5);

for j = 1 : 5
    
    y0 = [linspace(1e-3, -0.25, 10000), linspace(0.25, 4, 10000)];
    
    d = dlmread(filenameobs{j}); d = d(2);
    xMin = (-hp * y0 - (d - y0/2)) / sigmaErr(j);
    xMax = ( hp * y0 - (d - y0/2)) / sigmaErr(j);
        
    truth = exp(-y0.^2 / 2) .* (1 ./ y0) .* (normcdf(xMax) - normcdf(xMin));
    normPdf = truth / trapz(y0, truth);
 
    p(j) = plot(y0, truth / trapz(y0, truth), 'color', col(j, :));
    
%     fLim = exp(-y0.^2 / 2) .* (1 ./ y0) .* (y0 > 2 * sqrt(2) * d / (sqrt(2) + 1));
%     plot(y0, fLim ./ trapz(y0, fLim), ':', 'color', col(j, :));
    
    sigmaErrProb = sqrt(sigmaErr(j)^2  + 1 / 8);
    
%     truthMean = d / (2 * (sigmaErrProb^2 + 0.25));
%     truthStd = sigmaErrProb / sqrt(sigmaErrProb^2 + 0.25);
%     
%     exactMean = exp(-1/2) * d / (sigmaErr(j)^2 + exp(-1));
%     exactStd = sqrt(sigmaErr(j)^2 / (sigmaErr(j)^2 + exp(-1)));
    
%     y0 = linspace(-1, 4, 1000);
%     plot(y0, normpdf(y0, truthMean, truthStd), '--', 'color', col(j, :), 'linewidth', 1);

    meanex(j) = trapz(y0, y0 .* normPdf);
    secmomex(j) = trapz(y0, y0.^2 .* normPdf);
    varex(j) = secmomex(j) - meanex(j)^2;
    
    y_min = 2 * sqrt(2) * d / (sqrt(2) + 1);
    y_max = 2 * sqrt(2) * d / (sqrt(2) - 1);
    C = 0.5 * (expint(y_min^2 / 2) - expint(y_max^2 / 2));
    meanth(j) = sqrt(2 * pi) * (normcdf(y_max) - normcdf(y_min)) / C;
    secMomth(j) = (exp(-y_min^2/2) - exp(-y_max^2/2)) / C;
    varth(j) = secMomth(j) - meanth(j)^2;
        
end

dLim = exp(-1/2);
yMinLim = 2 * sqrt(2) * d / (sqrt(2) + 1);
yMaxLim = 2 * sqrt(2) * d / (sqrt(2) - 1);
CLim = 0.5 * (expint(yMinLim^2 / 2) - expint(yMaxLim^2 / 2));
meanthLim = sqrt(2 * pi) * (normcdf(yMaxLim) - normcdf(yMinLim)) / C;
secMomthLim = (exp(-yMinLim^2/2) - exp(-yMaxLim^2/2)) / C;
varthLim = secMomth(j) - meanth(j)^2;

p(j+1) = plot([1, 1], [0, max(normPdf)], 'k--', 'linewidth', 2);

for j = 1 : 5 
    plot([meanex(j), meanex(j)], [0, max(normPdf)], '--', 'color', col(j, :));
end

legend(p, {'$\sigma = 0.1$', '$\sigma = 0.05$', '$\sigma = 0.025$','$\sigma = 0.0125$', '$\sigma = 10^{-4}$', 'truth'}, 'interpreter', 'laTeX', 'location', 'NE');

yLim = 2 * sqrt(2) * d / (sqrt(2) + 1);
