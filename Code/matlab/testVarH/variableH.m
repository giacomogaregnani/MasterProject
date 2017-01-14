clc; clear; close all;

% ODE function f and initial condition
sigma = 10; rho = 28; beta = 8 / 3;
f = @(x) [sigma * (x(2) - x(1)); ...
    x(1) * (rho - x(3)) - x(2); ...
    x(1) * x(2) - beta * x(3)];
initialCond = [-10; -1; 40];
sizeSyst = 3;

% Final time
T = 100;

% Colours
lightgray = [0.7, 0.7, 0.7];

% Number of trajectories
nTrajectories = 5e3;

% Choose mean time step
hMean = 1e-2;

%% Random time stepping

expODE = zeros(nTrajectories, 1);
solutionODE = zeros(sizeSyst, nTrajectories);

parfor i = 1 : nTrajectories
    
    % Impose IC
    solution = initialCond;
    time = 0;
    
    while time(end) < T
        H = 2 * hMean * rand(1);
        newSol = solution(:, end) + H * f(solution(:, end));
        solution = [solution newSol];
        time = [time, time(end) + H];
    end
    
    if time(end) > T
        H = T - time(end - 1);
        newSol = solution(:, end-1) + H * f(solution(:, end-1));
        solution = [solution(:, 1:end-1) newSol];
        time = [time(1:end-1), time(end-1) + H];
    end
    
    %     plot(time, solution(1, :), 'color', lightgray);
    expODE(i) = solution(:, end)' * solution(:, end);
    solutionODE(:, i) = solution(:, end);
end

mExpODE = mean(expODE);

%% Associated SDE - Euler-Maruyama

expSDE = zeros(nTrajectories, 1);
solutionSDE = zeros(sizeSyst, nTrajectories);

parfor i = 1 : nTrajectories
    
    % Impose IC
    solution = initialCond;
    time = 0;
    
    while time(end) < T
        dW = sqrt(hMean) * randn(1);
        newSol = solution(:, end) + (hMean + sqrt(hMean/3) * dW) * f(solution(:, end));
        solution = [solution newSol];
        time = [time, time(end)+hMean];
    end
    
    %     plot(time, solution(1, :), 'color', lightgray);
    expSDE(i) = solution(:, end)' * solution(:, end);
    solutionSDE(:, i) = solution(:, end);
end
% ylim([-25, 25])

mExpSDE = mean(expSDE);

%% Probabilistic solver (Stuart Zyg.)

sigma = 0.5;
hMeanOneFive = hMean^(1.5);
noiseMultiplier = sqrt(sigma) * hMeanOneFive;

expPROB = zeros(nTrajectories, 1);
solutionPROB = zeros(sizeSyst, nTrajectories);

parfor i = 1 : nTrajectories
    
    % Impose IC
    solution = initialCond;
    time = 0;
    
    while time(end) < T
        xi = randn(sizeSyst, 1);
        newSol = solution(:, end) + hMean * f(solution(:, end)) + ...
            noiseMultiplier * xi;
        solution = [solution newSol];
        time = [time, time(end) + hMean];
    end
    
    expPROB(i) = solution(:, end)' * solution(:, end);
    solutionPROB(:, i) = solution(:, end);
end


%% Histogram of results at final time

figure
hold on
histogram(expODE, 'normalization', 'pdf', 'facealpha' , .5);
[f, xi] = ksdensity(expODE);
plot(xi, f)
xlim([0, 3000])
ylim([0, 1.2*1e-3])
figure
hold on
histogram(expSDE, 'normalization', 'pdf', 'facealpha' , .5);
[f, xi] = ksdensity(expSDE);
xlim([0, 3000])
ylim([0, 1.2*1e-3])
plot(xi, f)
figure
hold on
histogram(expPROB, 'normalization', 'pdf', 'facealpha' , .5);
[f, xi] = ksdensity(expPROB);
plot(xi, f)
xlim([0, 3000])
ylim([0, 1.2*1e-3])


%% Scatter plots
index = 1 : 5 : nTrajectories;
figure
scatter3(solutionODE(1, index), solutionODE(2, index), solutionODE(3, index), '.')
hold on
scatter3(solutionSDE(1, index), solutionSDE(2, index), solutionSDE(3, index), '.')
scatter3(solutionPROB(1, index), solutionPROB(2, index), solutionPROB(3, index), '.')
% xlim([-25, 25])







