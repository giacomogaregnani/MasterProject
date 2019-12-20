clc; clear; close all;

bound = 4;
pgon = polyshape([-bound -bound bound bound], [bound -bound -bound bound]);
tr = triangulation(pgon);
model = createpde();
tnodes = tr.Points';
telements = tr.ConnectivityList';
geometryFromMesh(model,tnodes,telements);
generateMesh(model, 'Hmax', 0.07);
pdemesh(model)
applyBoundaryCondition(model,'neumann','Edge',1:model.Geometry.NumEdges,'g',0);


p = @(x) sin(x);
dxp = @(x) cos(x);
dxxp = @(x) -sin(x);
eps = 0.1;
alpha = 1;
sigma = 0.5;

fcoef = @(location, state) -(alpha*location.x + 1/eps * dxp(location.x/eps)) .* state.ux ...
                           -(location.y - location.x) .* state.uy;
c = [sigma; 0];
a = @(location, state) alpha + 1/eps * dxxp(location.x/eps) + 1;
specifyCoefficients(model,'m',0,'d',0,'c',c,'a',a,'f',fcoef);

results=solvepde(model);
u = results.NodalSolution;
figure
pdeplot(model,'XYData',u)
title('Numerical Solution');
xlabel('x')
ylabel('y')