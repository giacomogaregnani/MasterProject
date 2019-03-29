function PlotVelocityField(Ux,Uy,delta)
% Plot of the PIECEWISE CONTINUOUS velocity field

xPlot = -1+delta/100:delta/10:1;
yPlot = -1+delta/100:delta/10:1;
[XX,YY] = meshgrid(xPlot,yPlot);

UxFunc = @(x,y) Ux(ceil((x+1)/delta),ceil((y+1)/delta));
UyFunc = @(x,y) Uy(ceil((x+1)/delta),ceil((y+1)/delta));

figure
surf(XX,YY,UxFunc(xPlot,yPlot)','EdgeColor','none')
figure
surf(XX,YY,UyFunc(xPlot,yPlot)','EdgeColor','none')

