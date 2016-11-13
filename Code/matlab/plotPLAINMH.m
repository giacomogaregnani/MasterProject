% read results
X = dlmread('../data/resultsMH.txt');
XRAM = dlmread('../data/resultsRAM.txt');

% contour plot of real distribution
f = @(x, y) exp(-10 * (x.^2 - y).^2 - (y - 0.25).^4);
xTest = -2 : 0.01 : 2;
yTest = -1 : 0.01 : 2;
[XX, YY] = meshgrid(xTest, yTest);

AX = [-2 2 -1.5 2];

contour(XX, YY, f(XX, YY), 'LineWidth', 2);
hold on
plot(X(:, 1), X(:, 2), 'r.')
axis(AX)

figure
contour(XX, YY, f(XX, YY), 'LineWidth', 2);
hold on
plot(XRAM(:, 1), XRAM(:, 2), 'r.')
axis(AX)
