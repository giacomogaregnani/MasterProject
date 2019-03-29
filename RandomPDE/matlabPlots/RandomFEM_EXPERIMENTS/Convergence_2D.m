clc; clear; close all;
%% Reference solution

uEx = [];
param = [];
zeroFunc = @(x,y) 0.*x.*y;

data   = read_DataFile('data_2D', 2, param);
data.diffusion = @(x, y, t, param) (x.*y).^0;
% data.diffusion = @(x,y,t,param)(1.3 + ... 1.2 * ((x - 0.25).^2 + (y - 0.25).^2 < 0.025) ...
%                                     - 0.5 * ((x - 0.75).^2 + (y - 0.75).^2 < 0.025));
data.diffusion = @(x,y,t,param)(0.5 * (x > 0.5) + 1.5 * (x < 0.5));                               
                             
data.param = param;
data.bcDir = @(x,y,t,param)(0*x.*y);
data.force = @(x,y,t,param) sin(2*pi*x) .* sin(2*pi*y);

% a = 5;
% data.force = @(x,y,t,param) -(a^2*exp(-a*((x - 1/2).^2 + (y - 1/2).^2)).*(2*y - 1).^2 - 2*a*exp(-a*((x - 1/2).^2 + (y - 1/2).^2)) ...
%     +a^2*exp(-a*((x - 1/2).^2 + (y - 1/2).^2)).*(2*x - 1).^2 - 2*a*exp(-a*((x - 1/2).^2 + (y - 1/2).^2)));
% p = 14;
% data.force = @(x,y,t,param)-(p^2*exp(p*x).*sin(2*pi*x).*sin(2*pi*y) + 4*p*pi*exp(p*x).*cos(2*pi*x).*sin(2*pi*y) - 4*pi^2*exp(p*x).*sin(2*pi*x).*sin(2*pi*y) ...
%     -4*pi^2*exp(p*x).*sin(2*pi*x).*sin(2*pi*y));

errData = data;
errData.uexact  = @(x,y,t,param) 0.*x.*y;
errData.uxexact = @(x,y,t,param) 0.*x.*y;
errData.uyexact = @(x,y,t,param) 0.*x.*y;

if isempty(uEx)
    % Space discretization
    [vertices, boundaries, elements] = poimesh('squareg', 300, 300);
    vertices = 0.5 * vertices + 0.5;
    
    [uEx, FESpaceRef, meshRef, ~, ~, solNorm]  = Elliptic_Solver(2, elements, vertices, boundaries, 'P1', errData);
    % Plot true solution
    FEMplot2D(meshRef, uEx, 'on', 45)
    title('true solution')
    
    % Continuous
    uExInt = pdeInterpolant(meshRef.vertices, meshRef.elements, uEx);
    uEx = @(x, y) uExInt.evaluate(x,y)';
    
    zeroFunc = @(x, y) 0.*x.*y;
else
    solNorm = integral2(uEx, 0, 1, 0, 1);
end


%% Convergence analysis

% Geometry specification [0,1]^2
g = [2 2 2 2
    0 1 1 0
    1 1 0 0
    1 1 0 0
    1 0 0 1
    0 0 0 0
    1 1 1 1];

pVec = 1.1:0.2:1.5;
NVec = 2.^[2:5];
hVec = 1 ./ NVec;
M = 6;
nP = length(pVec); nN = length(NVec);
errProbStruct = cell(nP, 1);
errInterpStruct = errProbStruct;
errProbH1Struct = errProbStruct;
errInterpH1Struct = errProbStruct;

col = colormap(parula(2*length(pVec)));
close
col = col(1:2:end, :);

for k = 1 : nP
    
    errL2 = zeros(M, length(hVec));
    errInterpL2 = errL2;
    errH1 = errL2;
    errInterpH1 = errL2;
    
    p = pVec(k);
    
    for i = 1 : nN
        
        h = hVec(i);
        [vert, bound, el] = initmesh(g,'hMax',h);
        
        [u, ~, meshDet, ~]  = Elliptic_Solver(2, el, vert, bound, 'P1', data);
        uFuncDet = pdeInterpolant(meshDet.vertices, meshDet.elements, u);
        data.uexact = @(x,y,t,param) evaluateFEMforError(uFuncDet,x,y,t,param);
        
        idxInternal = setdiff(1:size(meshDet.vertices, 2), [meshDet.boundaries(1, :), meshDet.boundaries(2, :)]);
        NVert = length(idxInternal);
        hPerturbation = zeros(1, NVert);
        itPert = 1;
        for ii = idxInternal
            [~, Neighbours] = find(meshDet.elements(1:3, :) == ii);
            hPerturbation(itPert) = min(meshDet.h(Neighbours));
            itPert = itPert + 1;
        end
        
        for jj = 1 : M
            meshProb = meshDet;
            meshProb.vertices(:, idxInternal) = meshProb.vertices(:, idxInternal) + hPerturbation.^p .* (rand(size(meshProb.vertices(:, idxInternal)))-0.5);
            [uProb, FESpaceProb, ~, ~, L2] = Elliptic_Solver(2, meshProb.elements, meshProb.vertices, meshProb.boundaries, 'P1', data, [], [], false);
            uFuncProb = pdeInterpolant(meshProb.vertices, meshProb.elements, uProb);
            errL2(jj, i) = L2;
            
            uInterp = uFuncDet.evaluate(meshProb.vertices(1, :), meshProb.vertices(2, :));
            uFuncInterp = pdeInterpolant(meshProb.vertices, meshProb.elements, uInterp);
            errInterpL2(jj, i) = FEM_error(uInterp,meshProb,data,FESpaceProb,[]);
            
            % Compute H1 Error
            uErr = uFuncDet.evaluate(meshRef.vertices(1, :), meshRef.vertices(2, :)) - uFuncProb.evaluate(meshRef.vertices(1, :), meshRef.vertices(2, :));
            [~, H1] = FEM_error(uErr, meshRef, errData, FESpaceRef, []);
            errH1(jj, i) = H1;
            uErrInterp = uFuncDet.evaluate(meshRef.vertices(1, :), meshRef.vertices(2, :)) - uFuncInterp.evaluate(meshRef.vertices(1, :), meshRef.vertices(2, :));
            [~, H1Interp] = FEM_error(uErrInterp, meshRef, errData, FESpaceRef, []);
            errInterpH1(jj, i) = H1Interp;
        end
    end
    
    errProbStruct{k} = mean(errL2);
    errInterpStruct{k} = mean(errInterpL2);
    errProbH1Struct{k} = mean(errH1);
    errInterpH1Struct{k} = mean(errInterpH1);
    
end

%%
figure
for k = 1 : length(pVec)
    loglog(hVec, errProbStruct{k}, 'o-', 'color', col(k, :))
    hold on
%     loglog(hVec, errInterpStruct{k}, 'x-', 'color', col(k, :))
    loglog(hVec, 2e-2 * hVec.^(pVec(k)+1), '--', 'color', col(k, :))
end
title('L2')

figure
for k = 1 : length(pVec)
    loglog(hVec, errProbH1Struct{k}, 'o-', 'color', col(k, :))
    hold on
%     loglog(hVec, errInterpH1Struct{k}, 'x-', 'color', col(k, :))
    loglog(hVec, 5e-2*hVec.^((pVec(k))/2), '--', 'color', col(k, :))
end
title('H1')

