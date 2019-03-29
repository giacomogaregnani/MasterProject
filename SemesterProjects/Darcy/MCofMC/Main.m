% Estimate TAU with DOUBLE MC simulation

clear
close all
clc

% Domain and general parameters
Time = [0, 20];
sigma = 0.01;
initialCondition = [-0.8; 0.8];
bounds = [-1, 1; -1, 1];
boundCond = 1;

% Parameters of outer MC
nRealizations = 100;
deltaU = 2 ^ -4;
pInlet = 1;
plotFields = 'False';
LMax = 8;
nu = 0.5;
LC = 0.05;
sigmaA = 1;

% Parameters of inner MC
N = Time(2) / deltaU;
nTrajCEM = 1e4;

% Initialization
nSigma = length(sigma);
exitTime = zeros(1, nSigma);

for j = 1 : nSigma
    tic
    for i = 1 : nRealizations
        
        % Generate the random field
        A = realizationRF(LMax(end), LC, nu, sigmaA, 1);
        NGridA = sqrt(length(A));
        A = reshape(A, NGridA, NGridA);
        
        % Solve Darcy
        [Ux, Uy] = SolveDarcyFF(A, pInlet, plotFields, deltaU);
        
        % Generate a BM
        W = BrownianMotion2D(Time, N, nTrajCEM);
        
        % Compute the exit time expectation
        currTau = CEMDarcy(initialCondition, sigma(j), bounds, boundCond, W, Time, Ux, Uy, deltaU);
        exitTime(j) = exitTime(j) + currTau;
        
        display(['sigma iteration: ', num2str(j), ' iterations remaining: ', num2str(nRealizations - i)])
    end
    clear W Ux Uy
    exitTime(j) = exitTime(j) / nRealizations;
    Time(j) = toc;
end
