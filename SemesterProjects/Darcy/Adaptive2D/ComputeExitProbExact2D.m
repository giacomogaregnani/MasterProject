function PhiExact = ComputeExitProbExact2D(Bounds, BoundCond, sigma, X0, Time)
% Compute reference EXIT PROBABILITY with PDE TOOLBOX

% Create the PDE model
model = createpde();

% Create the geometry and the mesh
R1 = [3,4,Bounds(1,1),Bounds(1,2),Bounds(2,2),Bounds(2,1),Bounds(1,2),Bounds(1,2),Bounds(1,1),Bounds(1,1)]';
g = decsg(R1);
geometryFromEdges(model,g);
MSH = generateMesh(model,'HMax',0.1);

% Set the coefficients of the equation
Coeff = specifyCoefficients(model,'m',0,'d',-1,'c',-1/2 * sigma^2,'a',0,'f',0);

% Set initial and boundary conditions
setInitialConditions(model,0);
if BoundCond == 0
    applyBoundaryCondition(model,'edge',1:4,'r',1);
else
    applyBoundaryCondition(model,'edge',[2,4],'r',1);
    applyBoundaryCondition(model,'edge',[1,3],'g',0);
end

% Solve and extract
results = solvepde(model,[Time(1):0.01:Time(2)]);
% figure
% u = results.NodalSolution;
% pdeplot(model,'xydata',u(:,end),'zdata',u(:,end),'colorbar','off','colormap','default')
% xlabel('x')
% ylabel('y')
% zlabel('\Phi(t=0)')

FinalTime = length(results.SolutionTimes);
PhiExact = interpolateSolution(results,X0(1),X0(2),FinalTime);

end