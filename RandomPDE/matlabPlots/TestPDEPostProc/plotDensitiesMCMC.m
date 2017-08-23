function plotDensitiesMCMC(sample, trueVal, ax, truePlot)

evalPoints = linspace(min(sample), max(sample), 400);
[dens, xi] = ksdensity(sample, evalPoints);

if nargin == 1

    plot(xi, dens)

elseif nargin == 2
    
    plot(xi, dens)
    hold on
    plot([trueVal, trueVal], [0, max(dens)], 'k--')
    
else
    
    if truePlot
        
        plot(ax, xi, dens)
        plot(ax, [trueVal, trueVal], [0, max(dens)], 'k--')
        
    else
        
        plot(ax, xi, dens)
        
    end
    
end