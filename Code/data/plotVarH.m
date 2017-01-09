% Plot some stuff

clc; clear; close all;

h = 0.025 * (0.9 .^ [0 : 86]);
hText = floor(h * 1e8);
nExp = length(h);

figure
hold on

% Read reference solution
yRefStruct = dlmread('fullResultsLorenz.txt');
yRef = yRefStruct(:, 2:end);
tRef = yRefStruct(:, 1);
clear yRefStruct;

plot(tRef, yRef(:, 1), 'k', 'LineWidth', 2);

lightgray = [0.6, 0.6, 0.6];

for i = 86 : 86
   
    y = dlmread(['fullTrajLorenz_h_', num2str(hText(i)), '.txt']);
    t = 0 : h(i) : 20 + h(i);
    
    plot(t, y(:, 1), 'color', lightgray);
    
end