function plotDensityIso(xPlot, yPlot, x, y, Dens)

[XX, YY] = meshgrid(xPlot, yPlot);

ZZ = griddata(x, y ,Dens, XX, YY);
ZZ(isnan(ZZ)) = 0;

contourf(XX, YY, ZZ);