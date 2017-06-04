clc; clear; close all;

%%
XDet = dlmread('invariantErrorDet.txt');
XProb = dlmread('invariantErrorProb.txt');

h = 0.1;
t = 0 : h : h * (size(XDet, 1) - 1);

%%
figure
plot(XProb(:, 3), XProb(:, 4));
xL = get(gca, 'xlim');
figure
plot(XDet(:, 3), XDet(:, 4));
set(gca, 'xlim', xL);

% Error on the Hamiltonian
figure
semilogy(t, abs(XDet(:, 5) - XDet(1, 5)));
hold on
semilogy(t, abs(XProb(:, 5) - XProb(1, 5)));
