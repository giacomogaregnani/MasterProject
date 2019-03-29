clc; clear; close all
%%

addpath('resultsOneStepVar');

sVec = []; hVec = [];

for i = 1 : 10
    figure
    hold on
    results = dlmread(['results', num2str(i-1), '.txt']); 
    
    t = results(:, 1); s = results(:, 2);
    sDelta = s(2:end) - s(1:end-1);
    
    plot(t(2:end), sDelta);
    
    sVec = [sVec, max(sDelta)];
    hVec = [hVec, t(2) - t(1)]; 
end

figure
loglog(hVec, sVec);
hold on
loglog(hVec, hVec.^(1.5))

order = log2(sVec(1:end-1) ./ sVec(2:end))