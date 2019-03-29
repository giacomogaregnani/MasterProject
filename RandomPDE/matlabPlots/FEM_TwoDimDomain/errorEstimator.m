clc; clear; close all;
%% Reference solution

% Space discretization
[vertices, boundaries, elements] = poimesh('squareg', 400, 400);
vertices = 0.5 * vertices + 0.5;

param = [];
data   = read_DataFile('InverseData', 2, param);
data.diffusion = @(x, y, t, param) (x.*y).^0;
% data.diffusion = @(x, y, t, param)(1.3 + 0.5 * ((x - 0.25).^2 + (y - 0.25).^2 < 0.025) ...
%                                   -0.5 * ((x - 0.75).^2 + (y - 0.75).^2 < 0.025));
% data.force = @(x,y,t,param)(8*pi^2*(sin(2*pi*x).*sin(2*pi*y)));
% data.force = @(x,y,t,param)(0*x.*y);
% data.force = @(x,y,t,param) exp(-1/1000 * ((x - 0.5).^2 + (y - 0.5).^2));
data.param = param;
data.bcDir = @(x,y,t,param)(0*x.*y);
data.bcNeu = @(x,y,t,param)(0*x.*y);
p = 18;
data.force = @(x,y,t,param) -(2^(4*p)*p*x.^p.*y.^p.*(p - 1).*(1 - x).^(p - 2).*(1 - y).^p ...
    - 2*2^(4*p)*p^2*x.^(p - 1).*y.^p.*(1 - x).^(p - 1).*(1 - y).^p + 2^(4*p)*p*x.^(p - 2).*y.^p.*(p - 1).*(1 - x).^p.*(1 - y).^p ... 
    + 2^(4*p)*p*x.^p.*y.^p.*(p - 1).*(1 - x).^p.*(1 - y).^(p - 2) - 2*2^(4*p)*p^2*x.^p.*y.^(p - 1).*(1 - x).^p.*(1 - y).^(p - 1) ...
    + 2^(4*p)*p*x.^p.*y.^(p - 2).*(p - 1).*(1 - x).^p.*(1 - y).^p);


% Solve
[uRef, ~, meshRef]  = Elliptic_Solver(2, elements, vertices, boundaries, 'P1', data);

% Plot true solution
FEMplot2D(meshRef, uRef)
title('true solution')

% Continuous
uRefInt = pdeInterpolant(meshRef.vertices, meshRef.elements, uRef);

%% Solution

M = 10;

for ii = 1 : 3
    N = 3 * 2^ii;
    h = 1 / N;
    [vert, bound, el] = poimesh('squareg', N, N);
    vert = 0.5 * vert + 0.5;
    
    % deterministic coarse solution
    [u, ~, meshCoarse, ~]  = Elliptic_Solver(2, el, vert, bound, 'P1', data);
    uInt = pdeInterpolant(meshCoarse.vertices, meshCoarse.elements, u);
    % deterministic error
%     errDet = abs(uInt.evaluate(meshRef.vertices(1, :), meshRef.vertices(2, :)) ...
%              -uRefInt.evaluate(meshRef.vertices(1, :), meshRef.vertices(2, :)));
    errDet = L2Error(meshCoarse, uRefInt, uInt);
%     plotP0(meshCoarse, L2Error(meshCoarse, uRefInt, uInt));
%     FEMplot2D(meshRef, err, 'off')

    
    % Perturb vertices
    errProb = 0;
    for jj = 1 : M
        MESH = meshCoarse;
        idxInternal = setdiff(1:size(MESH.vertices, 2), [MESH.boundaries(1, :), MESH.boundaries(2, :)]);
        MESH.vertices(:, idxInternal) = MESH.vertices(:, idxInternal) + 2 * h^2 * (rand(size(MESH.vertices(:, idxInternal)))-0.5);
        
        % Compute FEM solution
        U = Elliptic_Solver(2, MESH.elements, MESH.vertices, MESH.boundaries, 'P1', data, [], [], false);
        UInt = pdeInterpolant(MESH.vertices, MESH.elements, U);
        
        errProb = errProb + 1 / (M*h) * L2Error(meshCoarse, uInt, UInt);
%         errProb = errProb + 1 / (M*h) * abs(UInt.evaluate(meshRef.vertices(1, :), meshRef.vertices(2, :)) ...
%                                            -uInt.evaluate(meshRef.vertices(1, :), meshRef.vertices(2, :)));
    end
    plotP0(meshCoarse, errDet)
    zlim = caxis;
    plotP0(meshCoarse, errProb)
    caxis(zlim)
%     FEMplot2D(meshRef, errDet)
%     zlim = caxis;
%     FEMplot2D(meshRef, errProb)
%     caxis(zlim)
end
