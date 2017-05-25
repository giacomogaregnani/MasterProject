clc; clear; close all;

% ODE function f and initial condition
lambda = 1;
f = @(x) lambda * x;
initialCond = 1;
sizeSyst = 1;

% Final time
T = 1;

% Colours
lightgray = [0.7, 0.7, 0.7];

% Number of trajectories
nTrajectories = 1e4;

% Choose mean time step
hMean = 0.5;

% Solve for each time step with EE

expODE = zeros(nTrajectories, 1);
hVec = hMean;

for j = 1 : 7
    for i = 1 : nTrajectories
        
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
        
        expODE(i) = solution(end);
    end
    
    mExpODE(j) = mean(expODE);
    hMean = hMean / 2;
    hVec = [hVec hMean];
end
hVec = hVec(1:end-1);

%% SDE - analytical solution

expSDE = zeros(nTrajectories, 1);

hMean = 0.5;

for j = 1 : 7
    for i = 1 : 10^6
        
        W = sqrt(T) * randn(1);
        expSDE(i) = initialCond * exp((lambda - 0.5 * hMean / 3 * lambda^2) * T + ...
            sqrt(hMean/3) * lambda * W);
        
    end
    
    mExpSDE(j) = mean(expSDE);
    hMean = hMean / 2;
end

%% ERROR

err = abs(mExpODE - mExpSDE);
loglog(hVec, err, 'o-')
hold on
loglog(hVec, hVec, 'k--')

legend('weak error', 'slope 1', 'location', 'se')
xlabel('h')




