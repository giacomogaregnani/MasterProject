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
surf(xx, tt, y(:, 1:N)', 'edgecolor', [0.5, 0.5, 0.5])

