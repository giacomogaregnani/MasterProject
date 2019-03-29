clc; clear; close all
%% Reference solution

f = @(x) sin(2 * pi * x)';
xRef = linspace(0, 1, 10000);
uEx = @(x) 1 / (2 * pi)^2 * sin(2 * pi * x);

%% Convergence deterministic

N = 2.^[4:15];
hInt = 1 / 500;
xInt = linspace(0, 1, 501);
uExInt = uEx(xInt);

err = [];

for n = N
    h = 1 / n; 
    x = linspace(0, 1, n+1);
    A = 1 / h * gallery('tridiag', n-1);
    F = h * f(x(2:end-1));
    uVec = [0; (A \ F); 0];
    
    uInt = interp1(x, uVec, xInt);
    err = [err (trapz(xInt, (uInt - uExInt).^2))^0.5];
end

loglog(1./N, err, 'o-')
hold on
loglog(1./N, 1./(N.^2), 'k--')

%% Convergence probabilistic

errProb = [];
p = 1.2;

for n = N
    h = 1 / n; 
    x = linspace(0, 1, n+1);
    xProb = x;
    
    xProb(2:end-1) = x(2:end-1) + 0.5 * h^p * (-0.5 + rand(size(x(2:end-1))));
    
    A = assembleMatrixProb(xProb);
    F = assembleRHSProb(f, xProb);
    
    uVec = [0; (A \ F); 0];
    
    uInt = interp1(x, uVec, xInt);
    errProb = [errProb (trapz(xInt, (uInt - uExInt).^2))^0.5];
end

figure
loglog(1./N, errProb, 'o-')
hold on
loglog(1./N, 1./(N.^2), 'k--')

order = log2(errProb(1:end-1) ./ errProb(2:end))