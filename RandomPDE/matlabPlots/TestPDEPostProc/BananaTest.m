clc; clear; close all;
%%

B = 0.1;
d = 2;
nMCMC = 1e4;
alphaStar = 0.234;

if d > 2
    posterior = @(x) -x(1)^2 / 200 - 0.5 * (x(2) + B * x(1)^2 - 100 * B)^2 ...
        -0.5 * sum(x(3:end).^2);
elseif d == 2
    posterior = @(x) -x(1)^2 / 200 - 0.5 * (x(2) + B * x(1)^2 - 100 * B)^2;
end

[thetaAll, accRatio] = MCMCwithPost(posterior, ones(d, 1), nMCMC, alphaStar);

if d == 2
    dens = akde(thetaAll', thetaAll');
    
    xPlot = linspace(min(thetaAll(1, :)), max(thetaAll(1, :)), 1000);
    yPlot = linspace(min(thetaAll(2, :)), max(thetaAll(2, :)), 1000);
    [XX, YY] = meshgrid(xPlot, yPlot);
    ZZ = griddata(thetaAll(1, :), thetaAll(2, :), dens, XX, YY);
    ZZ(isnan(ZZ)) = 0;
    
    figure
    contourf(XX, YY, ZZ, 15)
    hold on
    scatter(thetaAll(1, :), thetaAll(2, :), 0.5, 'r.')
end