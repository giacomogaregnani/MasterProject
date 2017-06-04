%% Plot the stochastic process for modified equation
clc; clear; close all;

T = 5;
h = 0.01;
N = T / h;
avgPoly = 0;
avgTime = 0;
nMC = 10;
tPoly = linspace(0, T, T/h * 10);
T = linspace(0, N * h, N+1);

f1 = figure;
ax1 = axes('nextplot', 'add');
f2 = figure;
ax2 = axes('nextplot', 'add');

lightgrey = [0.6, 0.6, 0.6];

method = 'linear';

for j = 1 : nMC
    H = 2 * h * rand(1, N+1);
    
    if strcmp(method, 'poly')
        p = polyfit(T, H / h, N);
        thetaT = polyval(p, tPoly);
    elseif strcmp(method, 'linear')
        thetaT = interp1(T, H / h, tPoly, 'linear');
    end
    avgPoly = avgPoly + thetaT;
    
    plot(ax1, tPoly, thetaT, '--', 'color', lightgrey);
    plot(ax1, T, H/h, 'or');
    
    for i = 2 : length(thetaT)-1
        alphaT(i) = trapz(tPoly(1:i), thetaT(1:i));
    end
    plot(ax2, tPoly(2:end), alphaT, '--', 'color', lightgrey);
    avgTime = avgTime + alphaT;
    
end

avgPoly = avgPoly / nMC;
avgTime = avgTime / nMC;

plot(ax1, tPoly, avgPoly, 'b', 'linewidth', 2);
plot(ax1, [0, tPoly(end)], [1, 1], 'k', 'linewidth', 2);

plot(ax2, tPoly, tPoly, 'k', 'linewidth', 2);
plot(ax2, tPoly(2:end), avgTime, 'b', 'linewidth', 2);