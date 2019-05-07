clc; clear; close all
%%

addpath('RESULTSASSYR')

filename = {'MCMC_h_05_noise_1';
    'MCMC_h_05_noise_05';
    'MCMC_h_05_noise_025';
    'MCMC_h_05_noise_0125'};
filenameobs = {'observations1.txt', 'observations05.txt', ...
    'observations025.txt', 'observations0125.txt'};
enhanced = 1;
fontsizeLAB = getLatexTextSize('normalsize', 'enhanced', enhanced);
fontsizeTICK = getLatexTextSize('small', 'enhanced', enhanced);
W = 5; H = 5;

sigmaErr = [0.1, 0.05, 0.025, 0.0125];

col = colormap(lines(4));
close
col = col(1:end, :);

style = {':', '-.', '--', '-'};

%% ==========


fig = createFigure(W, H, 'enhanced',enhanced);
box on
hold on
for j = 1 : 4
    resultsProb = dlmread([filename{j}, '.txt']);
    [dens, xEvalProb] = ksdensity(resultsProb);
    
    %     p(j) = plot(xEval, dens, 'color', col(j, :));
    
    sigmaErrProb = sqrt(sigmaErr(j)^2  + 1 / 8);
    
    d = dlmread(filenameobs{j}); d = d(2);
    truthMean(j) = d / (2 * (sigmaErrProb^2 + 0.25));
    truthStd = sigmaErrProb / sqrt(sigmaErrProb^2 + 0.25);
    
    exactMean = exp(-1/2) * d / (sigmaErr(j)^2 + exp(-1));
    exactStd = sqrt(sigmaErr(j)^2 / (sigmaErr(j)^2 + exp(-1)));
    
    normPdfAN = normpdf(xEvalProb, truthMean(j), truthStd);
    p(j) = plot(xEvalProb, normPdfAN, style{j}, 'color', col(j, :));
    %     plot(xEval, normpdf(xEval, exactMean, exactStd), '--', 'color', col(j, :));
    
end

% for j = 1 : 4
%     plot([truthMean(j), truthMean(j)], [0, max(normPdfAN)], '--', 'color', col(j, :));
% end

p(j+1) = plot([1, 1], [0, max(dens)], 'k--');
% leg = legend([p(1), p(2), p(3), p(4)], {'$\sigma = 0.1$', '$\sigma = 0.05$', '$\sigma = 0.025$','$\sigma = 0.0125$'}, 'interpreter', 'laTeX', 'location', 'NW');
% columnlegend(2, {'$\sigma = 0.1$', '$\sigma = 0.05$', '$\sigma = 0.025$','$\sigma = 0.0125$'})

xlim([-1, 3])
xLim = get(gca, 'xLim');
xlabel('$\vartheta$', 'interpreter', 'latex')
set(gca, 'ytick', [])

title('Additive noise', 'interpreter', 'latex')

axpos = get(gca, 'innerposition');
print -depsc2 ../../../Reports/PaperRTSRK_18/VERSION15/ExPostAN.eps

%% ==========

fig = createFigure(W, H, 'enhanced',enhanced);
box on
hold on

hp = (0.5)^(1.5);
for j = 1 : 4
    
    resultsRTS = dlmread([filename{j}, '_rts.txt']);
    [dens, y0] = ksdensity(resultsRTS, 'numPoints', 200); %, 'bandwidth', 0.018);
    
    %     p(j) = plot(y0, dens, 'color', col(j, :));
    
    d = dlmread(filenameobs{j}); d = d(2);
    xMin = (-hp * y0 - (d - y0/2)) / sigmaErr(j);
    xMax = ( hp * y0 - (d - y0/2)) / sigmaErr(j);
    
    truth = exp(-y0.^2 / 2) .* (1 ./ y0) .* (normcdf(xMax) - normcdf(xMin));
    normPdfRTS = truth / trapz(y0, truth);
    p(j) = plot(y0, normPdfRTS, style{j}, 'color', col(j, :));
    
    mean(j) = trapz(y0, y0 .* normPdfRTS);
end

p(j+1) = plot([1, 1], [0, max(dens)], 'k--');
set(gca, 'xLim', xLim)

% for j = 1 : 4
%     plot([mean(j), mean(j)], [0, max(normPdfRTS)], '--', 'color', col(j, :));
% end

% leg = legend([p(1), p(2), p(3), p(4)], {'$\sigma = 0.1$', '$\sigma = 0.05$', '$\sigma = 0.025$','$\sigma = 0.0125$'}, 'interpreter', 'laTeX', 'location', 'NW');

xlabel('$\vartheta$', 'interpreter', 'latex')
set(gca, 'ytick', [])

title('Random time step', 'interpreter', 'latex')

set(gca, 'innerposition', axpos)
print -depsc2 ../../../Reports/PaperRTSRK_18/VERSION15/ExPostRTS.eps

%% ==========

fig = createFigure(W, H, 'enhanced',enhanced);
box on
hold on


for j = 1 : 4
    
    resultsDet = dlmread([filename{j}, '_det.txt']);
    [dens, y0] = ksdensity(resultsDet);
    
    %     p(j) = plot(xEval, dens, 'color', col(j, :));
    
    d = dlmread(filenameobs{j}); d = d(2);
    truthMean = d / (2 * (sigmaErr(j)^2 + 0.25));
    truthStd = sigmaErr(j) / sqrt(sigmaErr(j)^2 + 0.25);
    
    p(j) = plot(y0, normpdf(y0, truthMean, truthStd), style{j}, 'color', col(j, :));
    
end

p(j+1) = plot([1, 1], [0, max(dens)], 'k--');
set(gca, 'xLim', xLim)

% leg = legend([p(1), p(2), p(3), p(4)], {'$\sigma = 0.1$', '$\sigma = 0.05$', '$\sigma = 0.025$','$\sigma = 0.0125$'}, 'interpreter', 'laTeX', 'location', 'NW');

xlabel('$\vartheta$', 'interpreter', 'latex')
set(gca, 'ytick', [])

title('Deterministic', 'interpreter', 'latex')

set(gca, 'innerposition', axpos)
print -depsc2 ../../../Reports/PaperRTSRK_18/VERSION15/ExPostRK.eps

%% ==========

fig = createFigure(W, H, 'enhanced',enhanced);
box on
hold on

for j = 1 : 4
    
    y0 = linspace(0.5, 2, 1000);
    
    d = dlmread(filenameobs{j}); d = d(2);
    
    exactMean = exp(-1/2) * d / (sigmaErr(j)^2 + exp(-1));
    exactStd = sigmaErr(j) / sqrt(sigmaErr(j)^2 + exp(-1));
    
    p(j) = plot(y0, normpdf(y0, exactMean, exactStd), style{j}, 'color', col(j, :));
    
end

p(j+1) = plot([1, 1], [0, max(dens)], 'k--');
set(gca, 'xLim', xLim)

% leg = legend([p(1), p(2), p(3), p(4)], {'$\sigma = 0.1$', '$\sigma = 0.05$', '$\sigma = 0.025$','$\sigma = 0.0125$'}, 'interpreter', 'laTeX', 'location', 'NW');

xlabel('$\vartheta$', 'interpreter', 'latex')
set(gca, 'ytick', [])

title('True', 'interpreter', 'latex')

set(gca, 'innerposition', axpos)
print -depsc2 ../../../Reports/PaperRTSRK_18/VERSION15/ExPostTrue.eps

