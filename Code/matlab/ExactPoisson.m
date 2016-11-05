%% Load perturbed solution
clear; close all; clc

yStruct = dlmread('../data/refSolPoisson.txt');
T = yStruct(1, 1);
nData = yStruct(2, 1);
t = yStruct(3 : nData + 2, 1);
y = yStruct(nData + 3 : end, :);

%% MATRIX
A = -gallery('tridiag', nData);

%% Exact parameter value and distribution parameters
lambdaEx = 1; varLambda = 1; varData = 1e-3;

%% Build exact distribution

lambda = [0.5:0.001:1.5];
nLambda = length(lambda);
Q = ones(1, nLambda);
L = ones(1, nLambda);

for i = 1 : nLambda
    for j = 1 : T
        exSol = expm(lambda(i) * A * t(j)) * ones(nData, 1);
        % Prior (normal)
        Q(i) = Q(i) * 1 / (sqrt(2 * pi)) * exp(-(lambda(i) - lambdaEx)^2 / (2 * varLambda));
        % Likelihood
        L(i) = L(i) * (2 * pi * varData)^(-nData/2) * exp(-(exSol - y(j, :)')'*(exSol -  y(j, :)') / (2 * varData));
    end
end
% Exact normalized posterior distribution
P = Q .* L;
P = P / trapz(lambda, P);
varEx = var(lambda, P);

%% Load numerical solution Monte Carlo h variable N = 100

for i = 0 : 5
    resultsfile{i+1} = ['../data/testExact_MC_TenData_FixN_02_11_2016_01_28', ...
        num2str(i), '.txt'];
end
nExperience = length(resultsfile);

%% Plot results MonteCarlo h variable N = 100

% exact distribution
figure
hold on
plot(lambda, P, 'k', 'LineWidth', 2)

% numerical distribution h = 0.1
meanMCMC = zeros(nExperience, 1);
varMCMC = zeros(nExperience, 1);
for k = 1 : nExperience
    results = dlmread(resultsfile{k});
    x = results(10001:end);
    meanMCMC(k) = mean(x);
    varMCMC(k) = var(x);
    [f, xi] = ksdensity(x);
    plot(xi, f)
end
legend('Exact', 'h = 0.1', 'h = 0.05', 'h = 0.025', 'h = 0.0125', 'h = 0.0065', 'h = 0.0031')
set(gca, 'YTickLabel', '')
xlim([-1, 0])
title('Monte Carlo')

%% Load numerical solution GAUSS

for i = 0 : 4
    resultsfileGAUSS{i+1} = ['../data/testPoisson_GAUSS_03_11_2016_04_27', ...
        num2str(i), '.txt'];
end
nExperienceGAUSS = length(resultsfileGAUSS);

%% Plot results GAUSS

% exact distribution
figure
hold on
plot(lambda, P, 'k', 'LineWidth', 2)

% numerical distribution
meanGAUSS = zeros(nExperienceGAUSS, 1);
varGAUSS = zeros(nExperienceGAUSS, 1);
for k = 1 : nExperienceGAUSS
    results = dlmread(resultsfileGAUSS{k});
    x = results(10001:end);
    meanGAUSS(k) = mean(x);
    varGAUSS(k) = var(x);
    [f, xi] = ksdensity(x);
    plot(xi, f)
end
xlim([0.5 1.5])
legend('Exact', 'h = 0.1', 'h = 0.05', 'h = 0.025', 'h = 0.0125', 'h = 0.0065', 'h = 0.0031')
set(gca, 'YTickLabel', '')
title('Gauss')

%% Convergence WRT h
H_MC = zeros(nExperience, 1);
H_GAUSS = zeros(nExperienceGAUSS, 1);

for i = 1 : nExperience
    m = meanMCMC(i) - lambdaEx;
    V = varMCMC(i) + varEx;
    coeff = sqrt(2 * sqrt(varMCMC(i) * varEx) / V);
    H_MC(i) = 1 - coeff * exp(-0.25 * m^2 / V);
    H_MC(i) = sqrt(H_MC(i));
end

for i = 1 : nExperienceGAUSS
    m = meanGAUSS(i) - lambdaEx;
    V = varGAUSS(i) + varEx;
    coeff = sqrt(2 * sqrt(varGAUSS(i) * varEx) / V);
    H_GAUSS(i) = 1 - coeff * exp(-0.25 * m^2 / V);
    H_GAUSS(i) = sqrt(H_GAUSS(i));
end

h_TEST = 0.1 ./ 2.^[0:nExperience-1];
h_GAUSS = 0.1 ./ 2.^[0:nExperienceGAUSS-1];
h_ORDER = linspace(h_TEST(end), h_TEST(1), 1000);
figure
loglog(h_TEST, H_MC, 'o-')
hold on
loglog(h_GAUSS, H_GAUSS, 'o-')
loglog(h_ORDER, h_ORDER, 'k--')
legend('Monte Carlo', 'Gauss')
xlabel('h')
ylabel('err')

%% M2TIKZ OPTIONS
m2tikzoptions = ['ticklabel style={font=\Large},'...
    'legend style={font=\Large},' ...
    'title style={font=\Large}'];