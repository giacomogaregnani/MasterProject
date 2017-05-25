clear; clc; close all;

addpath('ResultsAdd', 'ResultsStep');

fileNameAddRoot = 'SolutionAdd';
fileNameStepRoot = 'SolutionStep';

info = dlmread('infoFile.txt');
nMC = info(1);
T = info(2);
h = info(3);
fullWrite = info(4);
timeVec = 0 : h : T;
finalTime = 200;
timeVec = timeVec(find(timeVec < finalTime + h));
endTime = length(timeVec);

% if variable fullWrite contained in the information file is true, then all
% the values of time integration are written in results files. Otherwise,
% only the solution at final time is recorded.


if fullWrite % All values, plot trajectories
    
    darkgrey = [0.5, 0.5, 0.5];
    
    enhanced = 1;
    W = 14;
    H = 5;
    
    fontsizelb = 'normalsize';
    fontsizetk = 'small';
    
    fig = createFigure(W, H, 'enhanced',enhanced);
    fontsizeLAB = getLatexTextSize(fontsizelb, 'enhanced', enhanced);
    fontsizeTICK = getLatexTextSize(fontsizetk, 'enhanced', enhanced);
    ax = axes();
        
    massStep = zeros(nMC, length(timeVec));
    
    for i = 0 : nMC - 1
        
        fileNameStep = [fileNameStepRoot, num2str(i), '.txt'];
        solStep = dlmread(fileNameStep);
        massStep(i + 1, :) = sum(solStep, 2);
        
        semilogy(timeVec, solStep(1:endTime, 3), 'color', darkgrey);
        hold on
        
    end
    
    ylim([1e-6, 1])
    yLimit = get(gca, 'yLim');
    
    ylabel('[X]', 'interpreter', 'LaTex');
    xlim([0, finalTime]);
    set(gca, 'xtick', [])
    set(gca, 'fontsize', fontsizeTICK);    
    set(get(gca, 'xlabel'), 'fontsize', fontsizeLAB);
    set(get(gca, 'ylabel'), 'fontsize', fontsizeLAB);
    
    print -depsc2 ../../Reports/RandTimeStep/VERSION3/OxPerox.eps
    
    fig = createFigure(W, H, 'enhanced',enhanced);
    
    massAdd = zeros(nMC, length(timeVec));
  
    for i = 0 : nMC - 1
        
        fileNameAdd = [fileNameAddRoot, num2str(i), '.txt'];
        solAdd = dlmread(fileNameAdd);
        massAdd(i + 1, :) = sum(solAdd, 2);

        
        semilogy(timeVec, solAdd(1:endTime, 3), 'color', darkgrey);
        hold on
        
    end
    
    xlabel('$t$', 'interpreter', 'LaTex');
    ylabel('[X]', 'interpreter', 'LaTex');
    xlim([0, finalTime]);
    
    set(gca, 'yLim', yLimit);

    set(gca, 'fontsize', fontsizeTICK);
    set(get(gca, 'xlabel'), 'fontsize', fontsizeLAB);
    set(get(gca, 'ylabel'), 'fontsize', fontsizeLAB);
    
    print -depsc2 ../../Reports/RandTimeStep/VERSION3/OxPeroxAdd.eps
    
    
else % Only final value, plot an histogram
    
    fileNameAdd = [fileNameAddRoot, '.txt'];
    solAddFull = dlmread(fileNameAdd);
    solAdd = sum(solAddFull.^2, 2);
    histogram(solAdd, 'normalization', 'pdf', 'facealpha' , .2)
    hold on
    
    
    fileNameStep = [fileNameStepRoot, '.txt'];
    solStepFull = dlmread(fileNameStep);
    solStep = sum(solStepFull.^2, 2);
    histogram(solStep, 'normalization', 'pdf', 'facealpha', .2)
    
end

