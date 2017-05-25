clc; clear; close all

addpath('ResultsFitznagBoxes')
nameRoot = 'testCubesD21H';
nameH = {'05', '1', '2', '3'};
result = zeros(length(nameH), 1); 

% compute volume of rectangle
x = dlmread([nameRoot, '3Fitznag.txt']);
vol = x(1, 3) * x(1, 4);

for i = 1 : length(nameH)
   [W, l] = computeDensity([nameRoot, nameH{i}, 'Fitznag'], 1e-1);
   result(i) = sum(vol * W);   
end