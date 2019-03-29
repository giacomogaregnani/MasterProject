function [Ux,Uy,dxA,time] = SolveDarcy(sigmaA,LC,nu,LMax,pInlet,plotfields)
% Solve UNCERTAIN DARCY's problem with PDE TOOLBOX (use Darcy.edp instead)

tic

% Generate the random field
A = realizationRF(LMax,LC,nu,sigmaA,1);
NGridA = sqrt(length(A));
A = reshape(A,NGridA,NGridA);
dxA = 2 / (NGridA - 1);
[X,Y] = meshgrid(-1:dxA:1,-1:dxA:1);

% Function random field (linear interpolation)
AFunc = @(region,state) interp2(X,Y,A,region.x,region.y);

% Create PDE model
model = createpde();
R1 = [3,4,-1,1,1,-1,1,1,-1,-1]';
g = decsg(R1);
geometryFromEdges(model,g);
MSH = generateMesh(model,'Hmax',0.03,'GeometricOrder','linear');

% Specify coefficients
Coeff = specifyCoefficients(model,'m',0,'d',0,'c',AFunc,'a',0,'f',0);
applyBoundaryCondition(model,'edge',[1,3],'g',0);
applyBoundaryCondition(model,'edge',4,'r',pInlet);
applyBoundaryCondition(model,'edge',2,'r',0);

results = solvepde(model);

if strcmp(plotfields,'True') == 1
    AInt = zeros(size(model.Mesh.Nodes,2),1);
    p = model.Mesh.Nodes;
    
    for i = 1 : length(AInt)
        AInt(i) = interp2(X,Y,A,p(1,i),p(2,i));
    end
    
    ux = -results.XGradients .* AInt;
    uy = -results.YGradients .* AInt;
    figure
    Pressure = results.NodalSolution;
    pdeplot(model,'xydata',Pressure,'zdata',Pressure,'colorbar','off','colormap','default')
    xlabel('x')
    ylabel('y')
    zlabel('p')
    
    figure
    pdeplot(model,'xydata',AInt,'zdata',AInt,'colorbar','off','colormap','default')
    xlabel('x')
    ylabel('y')
    zlabel('A')
    
    figure
    pdeplot(model,'xydata',ux,'zdata',ux,'colorbar','off','colormap','default')
    xlabel('x')
    ylabel('y')
    zlabel('u_x')
    
    figure
    pdeplot(model,'xydata',uy,'zdata',uy,'colorbar','off','colormap','default')
    xlabel('x')
    ylabel('y')
    zlabel('u_y')
end

% Definition of the results grid: use the A grid (with midpoints)
dxRes = dxA;
dyRes = dxA;
xRes = -1 + dxRes/2 : dxRes : 1 - dxRes / 2;
yRes = -1 + dyRes/2 : dyRes : 1 - dyRes / 2;

NxRes = length(xRes);
NyRes = length(yRes);

Axy = zeros(NxRes,NyRes);
Ux = Axy;
Uy = Axy;

for i = 1 : NxRes 
    for j = 1 : NyRes 
        Point = [-1+(i-0.5)*dxRes,-1+(j-0.5)*dyRes];
        Axy(i,j) = interp2(X,Y,A,Point(1),Point(2));
        [ux,uy] = evaluateGradient(results,Point(1),Point(2));
        Ux(i,j) = -Axy(i,j) * ux;
        Uy(i,j) = -Axy(i,j) * uy;
    end
end

time = toc;