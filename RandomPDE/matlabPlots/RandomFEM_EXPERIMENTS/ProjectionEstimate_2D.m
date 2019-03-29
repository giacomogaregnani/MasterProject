clc; clear; close all;
%% Reference solution

param = [];
data   = read_DataFile('data_2D', 2, param);
f = @(x,y) ( x .* (x < 0.5) + (1 - x) .* (x >= 0.5) ) .* sin(2*pi*y);

% Geometry specification [0,1]^2
g = [2 2 2 2
    0 1 1 0
    1 1 0 0
    1 1 0 0
    1 0 0 1
    0 0 0 0
    1 1 1 1];

pVec = 1:1;
NVec = 2.^[2:6];
hVec = 1 ./ NVec;
M = 1;
nP = length(pVec); 
nN = length(NVec);
errProb = cell(nP, 1);

for k = 1 : nP

    err = zeros(M, length(hVec));
    p = pVec(k);

    for i = 1 : nN
        
        h = hVec(i);
        
        [vert, bound, el] = initmesh(g,'hMax',h);
        mesh = buildMESH(2, el, vert, bound, 'P1', 4, data);
               
        F = f(mesh.vertices(1, :), mesh.vertices(2, :))';
        FInt = pdeInterpolant(mesh.vertices, mesh.elements, F);
        
        idxInternal = setdiff(1:size(mesh.vertices, 2), [mesh.boundaries(1, :), mesh.boundaries(2, :)]);
        NVert = length(idxInternal);
        hPerturbation = zeros(1, NVert);
        itPert = 1;
        for ii = idxInternal
            [~, Neighbours] = find(mesh.elements(1:3, :) == ii);
            hPerturbation(itPert) = min(mesh.h(Neighbours));
            itPert = itPert + 1;
        end
        
        for jj = 1 : M
            MESH = mesh;
            MESH.vertices(:, idxInternal) = MESH.vertices(:, idxInternal) + hPerturbation.^p .* (rand(size(MESH.vertices(:, idxInternal)))-0.5);
            FF = FInt.evaluate(MESH.vertices(1, :), MESH.vertices(2, :));
            FFInt = pdeInterpolant(MESH.vertices, MESH.elements, FF);
            err(jj, i) = sqrt(sum(L2Error(mesh, @(x,y)FInt.evaluate(x,y), @(x,y)FFInt.evaluate(x,y)).^2));
        end
        disp(['iteration ' num2str(k) '.' num2str(i) '/' num2str(length(pVec)) '.' num2str(length(hVec))])
        
    end
    
    if M > 1
        errProb{k} = mean(err);
    else
        errProb{k} = err;
    end
end

figure
for k = 1 : length(pVec)
    loglog(hVec, errProb{k}, 'o-')
    hold on
    loglog(hVec, hVec.^(pVec(k)+1), 'k--')
end
title('L2')

