%% Load perturbed solution
clear; close all; clc

yStruct = dlmread('../data/refSolTestOneDTenData.txt');
nData = yStruct(1);
T = yStruct(2);
t = yStruct(3 : nData + 2);
y = yStruct(nData + 3 : end);

%% Exact parameter value and distribution parameters
lambdaEx = -0.5; varLambda = 1; varData = 1e-3;

%% Build exact distribution

lambda = [-1 : 0.0001 : 0];
nLambda = length(lambda);
for i = 1 : nLambda
    exSol = exp(lambda(i) * t);
    % Prior (normal)
    Q(i) = 1 / (sqrt(2 * pi)) * exp(-(lambda(i) - lambdaEx)^2 / (2 * varLambda));
    % Likelihood
    L(i) = (2 * pi * varData)^(-nData/2) * exp(-sum((exSol - y).^2) / (2 * varData));
end
% Exact normalized posterior distribution
P = Q .* L;
P = P / trapz(lambda, P);
varEx = var(lambda, P);

%% Load numerical solution Monte Carlo h = 0.1, 0.05

for i = 0 : 4
    resultsfile{i+1} = ['../data/testExact_MC_01_11_2016_09_33', ...
        num2str(i), '.txt'];
    resultsfileh05{i+1} = ['../data/testExact_MC_h05_01_11_2016_10_18', ...
        num2str(i), '.txt'];
    resultsfileh01{i+1} = ['../data/testExact_MC_h01_01_11_2016_11_59', ...
        num2str(i), '.txt'];
    resultsfileh005{i+1} = ['../data/testExact_MC_h005_01_11_2016_05_08', ...
        num2str(i), '.txt'];
end
nExperience = length(resultsfile);

%% Plot results MonteCarlo h = 0.1, 0.05

% exact distribution
figure
hold on
plot(lambda, P, 'k', 'LineWidth', 2)

% numerical distribution h = 0.1
meanMCMC = zeros(nExperience, 1);
varMCMC = zeros(nExperience, 1);
col = {'b', 'r', 'm', 'c', 'g'};
for k = 1 : nExperience
    results = dlmread(resultsfile{k});
    x = results(10001:end);
    meanMCMC(k) = mean(x);
    varMCMC(k) = var(x);
    [f, xi] = ksdensity(x);
    plot(xi, f, col{k})
end
legend('Exact', 'N = 1', 'N = 10', 'N = 100', 'N = 1000', 'N = 10000')
set(gca, 'YTickLabel', '')
xlim([-1, 0])
title('h = 0.1')

% exact distribution
figure
hold on
plot(lambda, P, 'k', 'LineWidth', 2)

% numerical distribution h = 0.05
meanMCMCh05 = zeros(nExperience, 1);
varMCMCh05 = zeros(nExperience, 1);
col = {'b', 'r', 'm', 'c', 'g'};
for k = 1 : nExperience
    results = dlmread(resultsfileh05{k});
    x = results(10001:end);
    meanMCMCh05(k) = mean(x);
    varMCMCh05(k) = var(x);
    [f, xi] = ksdensity(x);
    plot(xi, f, col{k})
end
legend('Exact', 'N = 1', 'N = 10', 'N = 100', 'N = 1000', 'N = 10000')
set(gca, 'YTickLabel', '')
xlim([-1, 0])
title('h = 0.05')

% exact distribution
figure
hold on
plot(lambda, P, 'k', 'LineWidth', 2)

% numerical distribution h = 0.01
meanMCMCh01 = zeros(nExperience, 1);
varMCMCh01 = zeros(nExperience, 1);
col = {'b', 'r', 'm', 'c', 'g'};
for k = 1 : nExperience
    results = dlmread(resultsfileh01{k});
    x = results(10001:end);
    meanMCMCh01(k) = mean(x);
    varMCMCh01(k) = var(x);
    [f, xi] = ksdensity(x);
    plot(xi, f, col{k})
end
legend('Exact', 'N = 1', 'N = 10', 'N = 100', 'N = 1000', 'N = 10000')
set(gca, 'YTickLabel', '')
xlim([-1, 0])
title('h = 0.01')

% exact distribution
figure
hold on
plot(lambda, P, 'k', 'LineWidth', 2)

% numerical distribution h = 0.01
meanMCMCh005 = zeros(nExperience, 1);
varMCMCh005 = zeros(nExperience, 1);
col = {'b', 'r', 'm', 'c', 'g'};
for k = 1 : nExperience
    results = dlmread(resultsfileh005{k});
    x = results(10001:end);
    meanMCMCh005(k) = mean(x);
    varMCMCh005(k) = var(x);
    [f, xi] = ksdensity(x);
    plot(xi, f, col{k})
end
legend('Exact', 'N = 1', 'N = 10', 'N = 100', 'N = 1000', 'N = 10000')
set(gca, 'YTickLabel', '')
xlim([-1, 0])
title('h = 0.005')

%% Load numerical solution MLMC

for i = 0 : 2
    resultsfileMLMC{i+1} = ['../data/testML_Exact_02_11_2016_09_59', ...
        num2str(i), '.txt'];
end
nExperienceML = length(resultsfileMLMC);

%% Plot results MLMC

% exact distribution
figure
hold on
plot(lambda, P, 'k', 'LineWidth', 2)

% numerical distribution
meanMLMC = zeros(nExperienceML, 1);
varMLMC = zeros(nExperienceML, 1);
col = {'b', 'r', 'm', 'c', 'g'};
for k = 1 : nExperienceML
    results = dlmread(resultsfileMLMC{k});
    x = results(10001:end);
    meanMLMC(k) = mean(x);
    varMLMC(k) = var(x);
    [f, xi] = ksdensity(x);
    plot(xi, f, col{k})
end
xlim([-1 0])
H_ML = zeros(nExperienceML, 1);

for i = 1 : nExperienceML
    m = meanMLMC(i) - lambdaEx;
    V = varMLMC(i) + varEx;
    coeff = sqrt(2 * sqrt(varMLMC(i) * varEx) / V);
    H_ML(i) = 1 - coeff * exp(-0.25 * m^2 / V);
    H_ML(i) = sqrt(H_ML(i));
end

%% Load numerical solution GAUSS

for i = 0 : 2
    resultsfileGAUSS{i+1} = ['../data/testExact_GAUSS_TenData_02_11_2016_04_53', ...
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
xlim([-1 0])
H_GAUSS = zeros(nExperienceGAUSS, 1);

for i = 1 : nExperienceGAUSS
    m = meanGAUSS(i) - lambdaEx;
    V = varGAUSS(i) + varEx;
    coeff = sqrt(2 * sqrt(varGAUSS(i) * varEx) / V);
    H_GAUSS(i) = 1 - coeff * exp(-0.25 * m^2 / V);
    H_GAUSS(i) = sqrt(H_GAUSS(i));
end

%% Convergence WRT N
H_MC = zeros(nExperience, 1);
H_MC_h05 = H_MC;
H_MC_h01 = H_MC;
H_MC_h005 = H_MC;

for i = 1 : nExperience
    m = meanMCMC(i) - lambdaEx;
    V = varMCMC(i) + varEx;
    coeff = sqrt(2 * sqrt(varMCMC(i) * varEx) / V);
    H_MC(i) = 1 - coeff * exp(-0.25 * m^2 / V);
    H_MC(i) = sqrt(H_MC(i));
    m = meanMCMCh05(i) - lambdaEx;
    V = varMCMCh05(i) + varEx;
    coeff = sqrt(2 * sqrt(varMCMCh05(i) * varEx) / V);
    H_MC_h05(i) = 1 - coeff * exp(-0.25 * m^2 / V);
    H_MC_h05(i) = sqrt(H_MC_h05(i));
    m = meanMCMCh01(i) - lambdaEx;
    V = varMCMCh01(i) + varEx;
    coeff = sqrt(2 * sqrt(varMCMCh01(i) * varEx) / V);
    H_MC_h01(i) = 1 - coeff * exp(-0.25 * m^2 / V);
    H_MC_h01(i) = sqrt(H_MC_h01(i));
    m = meanMCMCh005(i) - lambdaEx;
    V = varMCMCh005(i) + varEx;
    coeff = sqrt(2 * sqrt(varMCMCh005(i) * varEx) / V);
    H_MC_h005(i) = 1 - coeff * exp(-0.25 * m^2 / V);
    H_MC_h005(i) = sqrt(H_MC_h005(i));
end

figure
N_THEORY = [1:1:10000];
N = 10.^[0:nExperience - 1];
loglog(N, H_MC, 'o-')
hold on
loglog(N, H_MC_h05, 'o-')
loglog(N, H_MC_h01, 'o-')
loglog(N, H_MC_h005, 'o-')
xlabel('N')
ylabel('error')
delta = 0.0;
loglog(N_THEORY, log(N_THEORY) ./ (N_THEORY.^(0.5 - delta)),'k--')
legend('h = 0.1', 'h = 0.05', 'h = 0.01', 'h = 0.005', 'theory')

%% Convergence WRT H
numberN = length(N);
numberH = 4;
H_MC_wrth = zeros(numberN, numberH);
for i = 1 : numberN
    H_MC_wrth(i, :) = [H_MC(i) H_MC_h05(i) H_MC_h01(i) H_MC_h005(i)];
end
h_TEST = [0.1, 0.05, 0.01, 0.005];
h_TEST_ML = [0.1, 0.05];
h_GAUSS = 0.1 ./ 2.^[0:nExperienceGAUSS-1];
figure
loglog(h_TEST, h_TEST, 'k--')
hold on
loglog(h_TEST, H_MC_wrth','o-')
loglog(h_TEST_ML, H_ML,'o-')
loglog(h_GAUSS, H_GAUSS, 'o-')
legend('theory', 'N = 1', 'N = 10', 'N = 100', 'N = 1000', 'N = 10000', 'MLMC')
xlabel('h')
ylabel('err')

%% M2TIKZ OPTIONS
m2tikzoptions = ['ticklabel style={font=\Large},'...
    'legend style={font=\Large},' ...
    'title style={font=\Large}'];