% Executable MAIN for one-dimensional case
clear
close all
clc

display('Running test for the 1D case')

% Define the problem
Time = [0, 5];
% SMOOTH f
V = @(x) 0.1 * (8 * x.^4 - 8 * x.^2 - x + 2);
dV = @(x) 0.1 * (32 * x.^3 - 16 * x - 1);
f = @(x) -dV(x);
g = @(x) 2;
X0 = 0;
Bounds = [-1, 1];
% CHANGE HERE TO [0, 0] FOR PURE KILLING BC's.
BoundCond = [0, 1];
N = 2.^[2 : 6];
M = 1e5;

% Compute the exact expectation of tau and the error
tauEx = ComputeExitTimeExact(X0,V,g,Bounds,BoundCond);
phiEx = ComputeExitProbFD(X0, Time, Bounds, BoundCond, f, g(1));

% Compute the BM
W = BrownianMotion(Time,N(end),M);

% Initialize
tauDEM = zeros(1,length(N));
tauCEM = tauDEM;
tNaive = tauDEM;
tBern = tauCEM;
phiDEM = tauDEM;
phiCEM = tauDEM;

for i = 1:length(N)
    % Compute the exit time and probability expectation
    [tauDEM(i),phiDEM(i),tNaive(i)] = ComputeExitTimeNaive(X0,f,g,Bounds,BoundCond,W(:,1:N(end)/N(i):end),Time);
    [tauCEM(i),phiCEM(i),tBern(i)] = ComputeExitTimeBernoulli(X0,f,g,Bounds,BoundCond,W(:,1:N(end)/N(i):end),Time);
    display([num2str(length(N) - i), ' iterations remaining'])
end

% Plot some trajectories
if BoundCond(2) == 0
    figure
    TrajectoriesForPlots(X0, f, g, Bounds, W(1:10,:), Time);
end

errDEMTau = abs(tauDEM - tauEx);
errCEMTau = abs(tauCEM - tauEx);
errDEMPhi = abs(phiDEM - phiEx);
errCEMPhi = abs(phiCEM - phiEx);

ErrorPlots(errCEMPhi, errDEMPhi, errCEMTau, errDEMTau, N, Time)

% Compute the orders
OrdersNaive = log2(errDEMTau(1:end-1)./errDEMTau(2:end));
OrdersBernoulli = log2(errCEMTau(1:end-1)./errCEMTau(2:end));
OrdersNaivePhi = log2(errDEMPhi(1:end-1)./errDEMPhi(2:end));
OrdersBernoulliPhi = log2(errCEMPhi(1:end-1)./errCEMPhi(2:end));

% Profiles of tau vs starting point
TauProfiles(V,dV,g,Bounds,BoundCond,W(1:1000,1:N(end)/N(end-1):end),Time,10);
PhiProfiles(f,g,Bounds,BoundCond,W(1:1000,1:N(end)/N(end-2):end),Time,10);

clear W
