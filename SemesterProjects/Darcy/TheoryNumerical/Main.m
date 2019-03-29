% Convergence wrt INTERPOLATION PARAMETERS 

clear
close all
clc

% Define the problem
Time = [0, 1];
sigma = 1;
g = @(x, y) sigma * eye(2);
f = @(x, y) 0.5 * [x^2 + y^2; abs(y - x)];
X0 = [0; 0];
Bounds = [-1, 1; -1, 1];
BoundCond = 1; % 0 for killing everywhere. 1 for two killing and two reflecting BCs. 2 for reflecting everywhere

% Define the size of the grid for "preprocessing" step
deltaU = 2 .^ [-1 : -1 : -8];

J = length(deltaU);
N_Balanced = Time(2) ./ deltaU;
N = N_Balanced(end) * 4;
M = 1e4;
I = length(N);

% Compute the BM
W = BrownianMotion2D(Time, N(end), M);

% Initialize
tauBernoulliC = zeros(1, length(N));
tBernoulliC = tauBernoulliC;
tauBernoulliL = zeros(1, length(N));
tBernoulliL = tauBernoulliL;
tauBernoulliCB = zeros(1, length(N));
tBernoulliCB = tauBernoulliCB;
XBernoulliC = zeros(2, length(N));
XBernoulliCB = zeros(2, length(N));
XBernoulliL = zeros(2, length(N));

for j = 1 : J
    
    [UxConst, UyConst] = constantInterpolation(f, deltaU(j));
    [tauBernoulliC(j), ~ ,tBernoulliC(j), XBernoulliC(:, j)] = ComputeExitTimeBernoulliConst(X0,g,Bounds,BoundCond,W(:,1 : N / ( 2*N_Balanced(end)) : end),Time,UxConst,UyConst,deltaU(j));
    [tauBernoulliCB(j), ~ ,tBernoulliCB(j), XBernoulliCB(:, j)] = ComputeExitTimeBernoulliConst(X0,g,Bounds,BoundCond,W(:,1:N / N_Balanced(j):end),Time,UxConst,UyConst,deltaU(j));
    clear UxConst UyConst
    
    %     [UxLin, UyLin] = linearInterpolation(f, deltaU(j));
    %     [tauBernoulliL(j), ~ ,tBernoulliL(j), XBernoulliL(j)] = ComputeExitTimeBernoulliLin(X0,g,Bounds,BoundCond,W(:,1 : N / N_Balanced(end) : end),Time,UxLin,UyLin,deltaU(j));
    %     clear UxLin UyLin
    
    display([num2str(J - j), ' iterations remaining'])
    
end

[RefTau, ~, tRef, XRef] = ComputeExitTimeBernoulli2D(X0, f, g, Bounds, BoundCond, W, Time);
errBernoulliTauC = abs(tauBernoulliC - RefTau);
errBernoulliTauCB = abs(tauBernoulliCB - RefTau);

for i = 1 : J
    errBernoulliXC(i) = norm(XBernoulliC(:, i) - XRef);
    errBernoulliXCB(i) = norm(XBernoulliCB(:, i) - XRef);
end

if BoundCond ~= 2
    % Error on tau
    figure
    loglog(deltaU, errBernoulliTauC, 'o-')
    hold on
    loglog(deltaU, errBernoulliTauCB, 'ro-')
    loglog(deltaU, deltaU, 'k--')
    grid on
    xlabel('\epsilon')
    legend('error', 'error_b', '\epsilon','Location','NW')
else
    % Error on the solution
    figure
    loglog(deltaU, errBernoulliXC, '*-')
    hold on
    loglog(deltaU, errBernoulliXCB, 'ro-')
    loglog(deltaU, deltaU * 4e-1, 'k--')
    grid on
    xlabel('\epsilon')
    legend('error', 'error_b', '\epsilon','Location','NW')
end

% ONLY IF LINEAR TESTS.
% figure
% loglog(deltaU, tBernoulliC, 'o-')
% hold on
% loglog(deltaU, tBernoulliCB, 'o-')
%
% figure
% loglog(deltaU, errBernoulliTauC, 'o-')
% hold on
% loglog(deltaU, errBernoulliTauL, 'o-')
% loglog(deltaU, deltaU, 'k--')
% loglog(deltaU, deltaU.^2, 'k')
% grid on
% xlabel('\epsilon')
% legend('errC', 'errL', '\epsilon', '\epsilon^2')
% figure
% loglog(deltaU, errBernoulliXC, 'o-')
% hold on
% loglog(deltaU, errBernoulliXCB, 'ro-')
% loglog(deltaU, errBernoulliXL, 'o-')
% loglog(deltaU, deltaU, 'k--')
% loglog(deltaU, deltaU.^2, 'k--')
% grid on
% xlabel('\epsilon')
% legend('error', 'error_b', 'error_l', '\epsilon', '\epsilon^2', 'Location','NW')
% title('error on X')

