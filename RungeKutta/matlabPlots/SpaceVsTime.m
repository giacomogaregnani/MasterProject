clc; clear; close all

addpath('ResultsSpaceVsTime');
 
% Ref solution
solTime = dlmread('SolutionTime.txt');

solSpaceRoot = 'SolutionSpace';
infoFileRoot = 'infoFile';

nExp = 3;

%%

for i = 1 : nExp
    
    figure
    hold on
    
    info = dlmread([infoFileRoot, num2str(i), '.txt']);
    T = info(2);
    h = info(3);
    timeVec = h : h : T;
    
    hText = floor(h * 1e4);
    
    solSpace = dlmread([solSpaceRoot, num2str(hText), '.txt']);
        
    plot(timeVec, solSpace)
    plot([h, T], [solTime(end), solTime(end)], 'k--', 'LineWidth', 2)
    
end