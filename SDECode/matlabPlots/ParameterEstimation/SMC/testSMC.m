clc; clear; close all;
%% Generate data

a = 1.0;
S = 0.5;
IC = 0;
T = 1;
N = 1000;

f = @(x) -a * x;
g = @(x) sqrt(2*S) * x.^0;

[t, obs] = EM_Traj(IC, f, g, T, N, 1);

% Perturb
sigmaNoise = 1e-3;
noisyObs = [IC, obs(2:end) + sigmaNoise * randn(1, N)];

%% SMC

aTest = linspace(0, 1, 100);
l = zeros(size(aTest));
M = 1000;


parfor i = 1 : length(aTest)    
    fTest = @(x) -aTest(i) * x;
    [~, ~, l(i)] = PF(IC, f, g, M, T, noisyObs, sigmaNoise);
    
%     % Choose the most likely path
%     idx = randsample(1:M, 1, true, w);
%     xChosen = x(idx, :);
%     
%     plot(t, xChosen)
%     hold on
%     plot(t, noisyObs)
    
end 

histogram(l)