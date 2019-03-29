clc; clear; close all;

%% Study the convergence of the probabilistic solver in an easy two-dimensional domain (with structured mesh)

%% REFERENCE SOLUTION

% Space discretization
[vertices, boundaries, elements] = poimesh('squareg', 500, 500);
nVertRef = size(vertices, 2);
hRef = 1 / 500;
vertices = 0.5 * vertices + 0.5;
param = [];
data   = read_DataFile('InverseData', 2, param);
data.diffusion = @(x, y, t, param) 1 + 0.*x.*y;
data.force = @(x,y,t,param)(8*pi^2*(sin(2*pi*x).*sin(2*pi*y)));
data.param = param;
data.bcDir = @(x,y,t,param)(0*x.*y);

% Solve
[uRef, ~, meshRef]  = Elliptic_Solver(2, elements, vertices, boundaries, 'P1', data);

% Plot true solution
FEMplot2D(meshRef, uRef)
title('true solution')

%% CONVERGENCE PROBABILISTIC

errProb = [];
nMC = 100;
p = 1.7;


for n = N
    
    h = 1 / n;
    display(['n = ' num2str(n)])
    
    [vert, bound, el] = poimesh('squareg', n, n);
    vert = 0.5 * vert + 0.5;
    [u, ~, mesh] = Elliptic_Solver(2, el, vert, bound, 'P1', data, [], [], false);

    MCMean = 0;
    parfor m = 1 : nMC
        MESH = mesh;
        idxInternal = setdiff(1:size(MESH.vertices, 2), [MESH.boundaries(1, :), MESH.boundaries(2, :)]);
        MESH.vertices(:, idxInternal) = MESH.vertices(:, idxInternal) ... 
                                        + 0.5 * h * (-0.5 + rand(size(MESH.vertices(:, idxInternal))));

        [uProb, ~, ~] = Elliptic_Solver(2, mesh.elements, MESH.vertices, mesh.boundaries, 'P1', data, [], [], false);
    
        uInt = pdeInterpolant(mesh.nodes, mesh.elements, uProb);
        uCoarseFine = evaluate(uInt, meshRef.nodes);
        
%         err = [err (1 / nVertRef * sum(uRef - uCoarseFine).^2)^0.5]; 
    
%         MCMean = MCMean + (1 / n^2 * sum(u - uProb).^2)^0.5;  

        MCMean = MCMean + (1 / nVertRef * sum(uRef - uCoarseFine).^2)^0.5;

    
    end
    
    errProb = [errProb MCMean/nMC];
    
end

figure
loglog(1./N, errProb, 'o-')
hold on
loglog(1./N, 1./(N.^p), 'k--')

order = log2(errProb(1:end-1) ./ errProb(2:end))