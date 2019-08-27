clc; clear; close all
%%
sol = dlmread('testParamSol.txt');
t = linspace(0, 100, size(sol, 2));
avg = dlmread('testParamTest.txt');

epsIdx = 4;
plot(t, sol(epsIdx, :))
hold on
plot(t, avg(epsIdx, :))

%%

names = {'testParam.txt', 'testParamAvg.txt'};
for i = 1 : 2    
    name = names{i};
    x = dlmread(name);
    
    hom = x(1, 1:2);
    x = x(2:end, :);
    n = size(x(:, 1));
    trueVals = [1, 0.5];
    
    figure
    hold on
    plot(x(:, 1), hom(1) * ones(n), 'k-.', 'linewidth', 1);
    plot(x(:, 1), trueVals(1) * ones(n), 'k--', 'linewidth', 1);
    for j = 2 : 2 : size(x, 2)
        plot(x(:, 1), x(:, j))
    end
    ylim = get(gca, 'ylim');
    set(gca, 'ylim', [0, ylim(2)])
    xlabel('$\epsilon$', 'interpreter', 'latex')
    ylabel('$A$', 'interpreter', 'latex')
    
    figure
    hold on
    plot(x(:, 1), hom(2) * ones(n), 'k-.', 'linewidth', 1);
    plot(x(:, 1), trueVals(2) * ones(n), 'k--', 'linewidth', 1);
    for j = 3 : 2 : size(x, 2)
        plot(x(:, 1), x(:, j))
    end
    ylim = get(gca, 'ylim');
    set(gca, 'ylim', [0, ylim(2)])
    xlabel('$\epsilon$', 'interpreter', 'latex')
    ylabel('$\Sigma$', 'interpreter', 'latex')
end