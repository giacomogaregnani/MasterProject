function tauExact = ComputeExitTimeExact2D(Bounds,BoundCond,sigma,X0)
% Compute EXIT TIME with PDE approach

model = createpde();
R1 = [3,4,Bounds(1,1),Bounds(1,2),Bounds(2,2),Bounds(2,1),Bounds(1,2),Bounds(1,2),Bounds(1,1),Bounds(1,1)]';
g = decsg(R1);
geometryFromEdges(model,g);
MSH = generateMesh(model,'Hmax',0.1);
Coeff = specifyCoefficients(model,'m',0,'d',0,'c',1,'a',0,'f',2/sigma^2);

if BoundCond == 0
    applyBoundaryCondition(model,'edge',1:4,'r',0);
    results = solvepde(model);
    tauExact = interpolateSolution(results,X0(1),X0(2));
else
    applyBoundaryCondition(model,'edge',[2,4],'r',0);
    applyBoundaryCondition(model,'edge',[1,3],'g',0);
    results = solvepde(model);
    tauExact = interpolateSolution(results,X0(1),X0(2));
end

figure
u = results.NodalSolution;
pdeplot(model,'xydata',u,'zdata',u,'colorbar','off','colormap','default')
xlabel('x')
ylabel('y')
zlabel('\tau')
title('Exit time on the whole domain - FEM')
end