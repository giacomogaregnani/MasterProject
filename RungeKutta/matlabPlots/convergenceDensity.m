clc; clear; close all

addpath('ResultsFitznagBoxes');
tol = 1e-1;

rootName = 'testCubesD21H';
hName = {'01', '05', '07', '08', '09', '1', '12', '15', '16', '2', '25', '3', '4'};

[wRef, l, xRef] = computeDensity([rootName, hName{1}, 'Fitznag'], tol);

err = zeros(length(hName) - 1, 1);
for i = 2 : length(hName)
   [W, l, x] = computeDensity([rootName, hName{i}, 'Fitznag'], tol);
   err(i - 1) = densDiff(xRef, wRef, x, W);
end

%% Plot
h = [0.05; 0.07; 0.08; 0.09; 0.1; 0.12; 0.15; 0.16; 0.2; 0.25; 0.3; 0.4];
loglog(h, err, 'o-')
hold on
loglog(h, h.^8, 'k--') 
