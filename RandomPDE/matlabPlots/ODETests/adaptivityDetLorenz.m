clc; clear; close all
%%

addpath('resultsAdaptDet')

nExp = 10;
filename = 'testAdapt';

figure
hold on

for i = 1 : nExp
    
    display(i)
    sol = dlmread([filename, num2str(i-1), '.txt']);
    plot(sol(:, 1), sol(:, 3), 'color', 0.6 * ones(3, 1))

end