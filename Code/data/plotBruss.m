% refSol
yStructRef = dlmread('refSolBruss.txt');
T = yStructRef(1, 1);
nData = yStructRef(2, 1);
t = yStructRef(3 : nData + 2, 1);
y = yStructRef(nData + 3 : end, :);

[~, N] = size(y);
N = N / 2;
x = linspace(0, 1, N);

[tt, xx] = meshgrid(t, x);
surf(tt, xx, y(:, 1:N)')

% numerical solution
yStructNum = dlmread('meanAndVariance.txt');
yNum = zeros(size(y));
var = yNum;
count = 1;
for i = 1 : 10
    yNum(i, :) = yStructNum(count, :)
    var(i, :) = diag(yStructNum(count + 1: count + 81 - 1, :));
    count = count + 81;
end

figure

surf(tt, xx, yNum(:, 1:N)')

firstSpec = yNum(:, 1);
firstVar = var(:, 1);
plot(t, firstSpec + firstVar);
hold on
plot(t, firstSpec + sqrt(firstVar));
plot(t, firstSpec - sqrt(firstVar));
