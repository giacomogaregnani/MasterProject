function p = plotDensitiesMCMC2(theta, col, ax)

t1m = min(theta(1, :)); t1M = max(theta(1, :));
t2m = min(theta(2, :)); t2M = max(theta(2, :));

xTest = linspace(t1m, t1M, 100);
yTest = linspace(t2m, t2M, 100);

[XX, YY] = meshgrid(xTest, yTest);

nKDE = floor(length(theta)^(1/2));

dens = akde(theta', [XX(:), YY(:)], nKDE);
ZZ = reshape(dens, 100, 100);

if nargin == 3
    p = contour(ax, XX, YY, ZZ, 10, 'color', col);
elseif nargin == 2
    p = contour(XX, YY, ZZ, 10, 'color', col);
elseif nargin == 1
    p = contour(XX, YY, ZZ, 10);
end

end