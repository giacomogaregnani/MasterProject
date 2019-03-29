clc; clear; close all
%%

x = dlmread('testParam.txt');

hom = x(1, 2:3);
x = x(2:end, :);
n = size(x(:, 1));

trueVals = [1, 0.5];

figure
hold on
plot(x(:, 1), hom(1) * ones(n), 'k-.', 'linewidth', 1);
plot(x(:, 1), trueVals(1) * ones(n), 'k--', 'linewidth', 1);
for j = 2 : 3 : size(x, 2)-2
   plot(x(:, 1), x(:, j)) 
end
ylim = get(gca, 'ylim');
set(gca, 'ylim', [0, ylim(2)])
xlabel('$\epsilon$', 'interpreter', 'latex')
ylabel('$A$', 'interpreter', 'latex')

figure
hold on
plot(x(:, 1), hom(1) * ones(n), 'k-.', 'linewidth', 1);
plot(x(:, 1), trueVals(1) * ones(n), 'k--', 'linewidth', 1);
for j = 4 : 3 : size(x, 2)
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
for j = 3 : 3 : size(x, 2)-1
   plot(x(:, 1), x(:, j)) 
end
ylim = get(gca, 'ylim');
set(gca, 'ylim', [0, ylim(2)])
xlabel('$\epsilon$', 'interpreter', 'latex')
ylabel('$\Sigma$', 'interpreter', 'latex')