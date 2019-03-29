% verify the convergence of DEM and CEM in a square domain. 

clear
close all
clc

display('Running test for the 2D case')

% Solution of the transport diffusion SDE with velocity field v = -Ax;

% Define the problem 
Time = [0,3];
sigma = 1;

V = @(x,y) zeros(2,1) * x * y;
dV = @(x,y) zeros(2,1) * x * y;
f = @(x,y) -dV(x,y);
g = @(x,y) sigma * eye(2);
X0 = [0;0];
Bounds = [-1,1;-1,1];
BoundCond = 0; % 0 for killing everywhere. 1 for two killing and two reflecting BCs.
N = 2.^[3:7];
M = 1e4;

% Compute the BM
W = BrownianMotion2D(Time,N(end),M);

% Initialize
tauNaive = zeros(1,length(N));
tauBernoulli = tauNaive;
phiNaive = tauNaive;
phiBernoulli = tauNaive;
tNaive = tauNaive;
tBernoulli = tauNaive;

for i = 1:length(N)
    % Compute the exit time expectation
    [tauNaive(i),phiNaive(i),tNaive(i)] = ComputeExitTimeNaive2D(X0,f,g,Bounds,BoundCond,W(:,1:N(end)/N(i):end),Time);
    [tauBernoulli(i),phiBernoulli(i),tBernoulli(i)] = ComputeExitTimeBernoulli2D(X0,f,g,Bounds,BoundCond,W(:,1:N(end)/N(i):end),Time);
    display([num2str(length(N) - i), ' iterations remaining'])
end

% Compute the exact expectation of tau and the error
tauEx = ComputeExitTimeExact2D(Bounds,BoundCond,sigma,X0);
phiEx = ComputeExitProbExact2D(Bounds,BoundCond,sigma,X0,Time);
errDEMTau = abs(tauNaive - tauEx);
errCEMTau = abs(tauBernoulli - tauEx);
errDEMPhi = abs(phiNaive - phiEx);
errCEMPhi = abs(phiBernoulli - phiEx);

% Plot the error for orders analysis on Tau
h = (Time(2)-Time(1))./N;
IndForPlots = ceil(length(N)/2);
figure
loglog(h,errDEMTau,'ro-')
hold on
loglog(h,errCEMTau,'b*-')
loglog(h,sqrt(h)*(errDEMTau(IndForPlots)/sqrt(h(IndForPlots))),'k--')
loglog(h,h*(errCEMTau(IndForPlots)/h(IndForPlots)),'k')
grid on
h_legend = legend('DEM','CEM','h^{0.5}','h');
set(h_legend,'Location','best','FontSize',13);
xlabel('h')
title('Convergence of the mean exit time')

% Plot the error for orders analysis on Phi
% figure
% loglog(h,errDEMPhi,'ro-')
% hold on
% loglog(h,errCEMPhi,'b*-')
% loglog(h,sqrt(h)*(errDEMPhi(IndForPlots)/sqrt(h(IndForPlots))),'k--')
% loglog(h,h*(errCEMPhi(IndForPlots)/h(IndForPlots)),'k')
% grid on
% h_legend = legend('err_h^{d,\Phi}','err_h^{c,\Phi}','h^{0.5}','h');
% set(h_legend,'Location','northwest','FontSize',13);
% xlabel('h')

% Compute the orders
OrdersNaiveTau = log2(errDEMTau(1:end-1)./errDEMTau(2:end));
OrdersBernoulliTau = log2(errCEMTau(1:end-1)./errCEMTau(2:end));
OrdersNaivePhi = log2(errDEMPhi(1:end-1)./errDEMPhi(2:end));
OrdersBernoulliPhi = log2(errCEMPhi(1:end-1)./errCEMPhi(2:end));

clear W

