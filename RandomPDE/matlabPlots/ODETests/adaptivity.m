clc; clear; close all
%%
addpath('resultsAdapt')

nSimPlot = 4;
sol = dlmread(['testAdapt', num2str(nSimPlot), '.txt']);
solEmb = dlmread(['testAdapt', num2str(nSimPlot), 'Emb.txt']);
err = dlmread('testAdaptsummary.txt');
errEmb = err(1:2:end, :);
errProb = err(2:2:end, :);

meanHEmb = errEmb(:, 5)'
meanHProb = errProb(:, 5)'
rejEmb = errEmb(:, 4)'
rejProb = errProb(:, 4)'
nFuncEvalEmb = errEmb(:, 3)'
nFuncEvalProb = errProb(:, 3)'
errorEmb = errEmb(:, 2)'
errorProb = errProb(:, 2)'
errorProb2 = errProb(:, end-1)'
stddevProb = errProb(:, end)'
tol = errProb(:, 1);

d = 2;
comp = 1;

subplot(2, 2, 1)
plot(sol(:, 1), sol(:, 1+comp:d:end-1), 'color', 0.6 * ones(3, 1))
hold on
plot(solEmb(:, 1), solEmb(:, 1+comp:d:end-1), 'k')

subplot(2, 2, 3)
semilogy(sol(2:end,1), sol(2:end, 1) - sol(1:end-1, 1), 'color', 0.6 * ones(3, 1))
hold on
semilogy([sol(2,1), sol(end,1)], [meanHProb(nSimPlot+1), meanHProb(nSimPlot+1)], 'color', 0.6 * ones(3, 1))
semilogy(solEmb(2:end,1), solEmb(2:end, 1) - solEmb(1:end-1, 1), 'k')
semilogy([solEmb(2,1), solEmb(end,1)], [meanHEmb(nSimPlot+1), meanHEmb(nSimPlot+1)], 'k')

subplot(2, 2, 2)
loglog(tol, errorEmb, 'k*-')
hold on
loglog(tol, errorProb, 'o-', 'color', 0.6 * ones(3, 1))
loglog(tol, errorProb2, 'o--', 'color', 0.6 * ones(3, 1))
loglog(tol, tol, 'k')
set(gca, 'xdir', 'reverse')
xlabel('tol')
ylabel('err')
legend('emb', 'prob')

subplot(2, 2, 4)
loglog(errorEmb, nFuncEvalEmb, 'k*-')
hold on
loglog(errorProb, nFuncEvalProb, 'o-', 'color', 0.6 * ones(3, 1))
loglog(errorProb2, nFuncEvalProb, 'o--', 'color', 0.6 * ones(3, 1))
set(gca, 'xdir', 'reverse')
xlabel('err')
ylabel('work')
legend('emb', 'prob')

figure
hold on
plot(sol(:, 1), sol(:, end), 'color', 0.6 * ones(3, 1)) 
plot(solEmb(:, 1), solEmb(:, end), 'k')
plot([sol(1, 1), sol(end, 1)], [tol(nSimPlot+1), tol(nSimPlot+1)], 'k--')