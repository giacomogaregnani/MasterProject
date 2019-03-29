clc; clear; close all

filename = 'calibration';

figure
hold on
for i = 1 : 4
    
    y = dlmread([filename, num2str(i), 'Growth2.txt']);
    N = size(y, 1);
    
    plot(y(1:end, 1), y(1:end, 2))
    
    variance(i) = y(end, 2);
    
end

figure
h = [0.2, 0.1, 0.05, 0.025];
loglog(h, variance, 'o-');
hold on
loglog(h, 10*h.^1.5, 'k--');