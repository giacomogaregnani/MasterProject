clc; clear; close all

%% H = 0.001

fileName = 'mcmcFITZNAG06_12_2016_03_58_';

results = dlmread([fileName, '.txt']);
N0001 = results(:, 1);
trueVal = [0.2, 0.2, 3.0] * [0.2, 0.2, 3.0]';
completeTheta0001 = dlmread([fileName, num2str(N0001(end)), '.txt']);
cumMean0001 = cumsum(completeTheta0001) ./ [1 : N0001(end)]';

cumVar001 = zeros(N0001, 1);
parfor i = 1 : N0001
   cumVar001(i) = var(completeTheta0001(1 : i)); 
end


%% H = 0.001

fileName = 'mcmcFITZNAG06_12_2016_02_58_';

results = dlmread([fileName, '.txt']);
N001 = results(:, 1);
trueVal = [0.2, 0.2, 3.0] * [0.2, 0.2, 3.0]';
completeTheta001 = dlmread([fileName, num2str(N001(end)), '.txt']);
cumMean001 = cumsum(completeTheta001) ./ [1 : N001(end)]';

cumVar001 = zeros(N001, 1);
parfor i = 1 : N001
   cumVar001(i) = var(completeTheta001(1 : i)); 
end


%% H = 0.01
fileName = 'mcmcFITZNAG06_12_2016_02_32_';

results = dlmread([fileName, '.txt']);
N01 = results(:, 1);

trueVal = [0.2, 0.2, 3.0] * [0.2, 0.2, 3.0]';

completeTheta = dlmread([fileName, num2str(N01(end)), '.txt']);
cumMean01 = cumsum(completeTheta) ./ [1 : N01(end)]';

cumVar01 = zeros(N001, 1);
parfor i = 1 : N001
   cumVar01(i) = var(completeTheta(1 : i)); 
end


%% H = 0.05
fileName = 'mcmcFITZNAG06_12_2016_03_31_';

results = dlmread([fileName, '.txt']);
N1 = results(:, 1);
trueVal = [0.2, 0.2, 3.0] * [0.2, 0.2, 3.0]';
completeTheta = dlmread([fileName, num2str(N1(end)), '.txt']);
cumMean05 = cumsum(completeTheta) ./ [1 : N1(end)]';

cumVar05 = zeros(N001, 1);
parfor i = 1 : N001
   cumVar05(i) = var(completeTheta(1 : i)); 
end

%% PLOOOOOT
figure
hold on
plot(completeTheta001, 'color', [0.65 0.65 0.65])
plot(cumMean001, 'g', 'LineWidth', 2)
plot(cumMean01, 'b', 'LineWidth', 2)
plot(cumMean05, 'r', 'LineWidth', 2)
plot([1, N01(end)], [trueVal, trueVal], 'k--', 'LineWidth', 2)
legend('g(\theta_k)', 'h = 0.001', 'h = 0.01', 'h = 0.05', 'g(\theta)')

figure
hold on
plot(cumVar001, 'g', 'LineWidth', 2)
plot(cumVar01, 'b', 'LineWidth', 2)
plot(cumVar05, 'r', 'LineWidth', 2)
legend('h = 0.001', 'h = 0.01', 'h = 0.05')