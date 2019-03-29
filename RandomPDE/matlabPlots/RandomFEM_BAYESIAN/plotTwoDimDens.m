function plotTwoDimDens(sample, nPoints, col)
% Sample must be a two columns vector

if nargin < 2
    nPoints = 50;
end

xEval = linspace(min(sample(:, 1))-0.1, max(sample(:, 1))+0.1, nPoints);
yEval = linspace(min(sample(:, 2))-0.1, max(sample(:, 2))+0.1, nPoints);
[XX, YY] = meshgrid(xEval, yEval);

[n, d] = size(sample);
bandwidth = 2 * std(sample) * ((4 / (n * (d + 2)))^(1 / (d + 4)));
m = mean(sample);

dens = mvksdensity(sample, [XX(:), YY(:)], 'Bandwidth', bandwidth);
dens = reshape(dens, nPoints, nPoints);

if nargin == 3
    contour(XX, YY, dens, 10, 'color', col)
else
    contour(XX, YY, dens, 10)
end

return