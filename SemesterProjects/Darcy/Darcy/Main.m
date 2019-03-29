% Executable MAIN for the DARCY PROBLEM and STOCHASTIC PARTICLE TRANSPORT

clear
close all
clc

% Define the problem
Time = [0,10];
sigma = 0.3;
g = @(x,y) sigma * eye(2);
X0 = [-0.8; -0.8];
Bounds = [-1,1;-1,1];
BoundCond = 1;

% Define the number of points for the grids for A
display('Generating random field')
LMax = 7;
nu = 0.5;
sigmaA = 1;
LC = 0.05;
A = realizationRF(LMax, LC, nu, sigmaA, 1);
NGridA = sqrt(length(A));
A = reshape(A, NGridA, NGridA);

% Darcy parameters
pInlet = 1;
plotfields = 'True';

% Define the size of the grid for "preprocessing" step
deltaU = 2^-6;
J = length(deltaU);
N = Time(2) / deltaU;
M = 1e3;

% Compute the BM
W = BrownianMotion2D(Time, N, M);

% Solve Darcy
display('Solving Darcy equation with FreeFem++')
[Ux,Uy] = SolveDarcyFF(A,pInlet,plotfields,deltaU);

% Compute the exit time expectation
display('Performing CEM')
[tauCEM, phiCEM, tCEM] = ComputeExitTimeBernoulliDarcy(X0, g, Bounds, BoundCond, W, Time, Ux, Uy, deltaU);
[ExpTau,ExpPhi,t] = ComputeExitTimeNaiveDarcyPlot(X0, g, Bounds, BoundCond, W(1:40,:), Time, Ux, Uy, deltaU);

display(['Mean exit time: ', num2str(tauCEM)])

clear W
