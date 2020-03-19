clc; clear; close all;
%% Plot parameters
enhanced = 1;
fontsizeLAB = getLatexTextSize('normalsize', 'enhanced', enhanced);
fontsizeTICK = getLatexTextSize('small', 'enhanced', enhanced);

W = 6.7; H = 6.7;

%% Generate reference mesh (long)

% Geometry specification [0,1]^2

% Unit square
g = [2 2 2 2
    0 1 1 0
    1 1 0 0
    1 1 0 0
    1 0 0 1
    0 0 0 0
    1 1 1 1];

% L-shaped domain
% g = [2   2   2  2  2   2
%      1/2 1   1  0  0   1/2
%      1   1   0  0  1/2 1/2
%      1/2 1/2 1  1  0   0
%      1/2 1   1  0  0   1/2
%      1   1   1  1  1   1
%      0   0   0  0  0   0  ];

[vertices, boundaries, elements] = initmesh(g, 'hMax', 1/100);

%% Reference solution


uEx = [];
param = [];
zeroFunc = @(x,y) 0.*x.*y;

data = read_DataFile('data_2D', 2, param);
% data.diffusion = @(x, y, t, param) (x.*y).^0;
% data.diffusion = @(x,y,t,param)(1 - 0.8 * ((x - 0.75).^2 + (y - 0.75).^2 < 0.04));
data.diffusion = @(x,y,t,param) (0.1 * (x > 0.5) + 0.8 * (x <= 0.5)).*(y.^0);

data.param = param;
data.bcDir = @(x,y,t,param)(0*x.*y);
data.bcNeu = @(x,y,t,param)(0*x.*y);
% data.bcNeu = @(x,y,t,param)(x).^2;

% a = 200;
% data.force = @(x,y,t,param) -(a^2*exp(-a*((x - 1/2).^2 + (y - 1/2).^2)).*(2*y - 1).^2 - 2*a*exp(-a*((x - 1/2).^2 + (y - 1/2).^2)) ...
%     +a^2*exp(-a*((x - 1/2).^2 + (y - 1/2).^2)).*(2*x - 1).^2 - 2*a*exp(-a*((x - 1/2).^2 + (y - 1/2).^2)));
a = 4;
data.force = @(x,y,t,param)(100 * (sin(a*pi*x).*sin(a*pi*y)));
% p = 42;
% data.force = @(x,y,t,param) -(2^(4*p)*p*x.^p.*y.^p.*(p - 1).*(1 - x).^(p - 2).*(1 - y).^p ...
%     - 2*2^(4*p)*p^2*x.^(p - 1).*y.^p.*(1 - x).^(p - 1).*(1 - y).^p + 2^(4*p)*p*x.^(p - 2).*y.^p.*(p - 1).*(1 - x).^p.*(1 - y).^p ...
%     + 2^(4*p)*p*x.^p.*y.^p.*(p - 1).*(1 - x).^p.*(1 - y).^(p - 2) - 2*2^(4*p)*p^2*x.^p.*y.^(p - 1).*(1 - x).^p.*(1 - y).^(p - 1) ...
%     + 2^(4*p)*p*x.^p.*y.^(p - 2).*(p - 1).*(1 - x).^p.*(1 - y).^p);
% p = 14;
% data.force = @(x,y,t,param)-(p^2*exp(p*x).*sin(2*pi*x).*sin(2*pi*y) + 4*p*pi*exp(p*x).*cos(2*pi*x).*sin(2*pi*y) - 4*pi^2*exp(p*x).*sin(2*pi*x).*sin(2*pi*y) ...
%                             -4*pi^2*exp(p*x).*sin(2*pi*x).*sin(2*pi*y));

% Space discretization
[uExNodes, ~, meshRef]  = Elliptic_Solver(2, elements, vertices, boundaries, 'P1', data);

% Plot true solution
fig = createFigure(W, H, 'enhanced',enhanced);
FEMplot2D(meshRef, uExNodes, 'on', 10, 201, 10)
colorbar off
title({'reference solution'}, 'interpreter', 'latex')
xlabel('$x_1$', 'interpreter', 'latex')
ylabel('$x_2$', 'interpreter', 'latex')
axis square
inPos = get(gca, 'innerpos');
% export_fig(fig, '../../../Reports/DraftPDE_18/VERSION7/Figures/MeshAdapt2D_ReferenceSol.eps', '-nocrop', '-painters')

% Continuous
uExInt = pdeInterpolant(meshRef.vertices, meshRef.elements, uExNodes);
uEx = @(x, y) uExInt.evaluate(x,y);
data.uexact = @(x,y,t,param) 0.*x.*y;
[~, ~, ~, ~, solNorm, ~, ~] = Elliptic_Solver(2, elements, vertices, boundaries, 'P1', data);
data.uexact = @(x,y,t,param) evaluateFEMforError(uExInt,x,y,t,param);


%% Solution

N = 5;
[vert, bound, el] = initmesh(g, 'hMax', 1 / N);

M = 10;
flag = true;

tol = 1e-1;
iter = 0;
errEst = []; errTrue = []; NVec = [];

while flag
    
    iter = iter + 1;
    
    % deterministic coarse solution
    data.uexact = @(x,y,t,param) evaluateFEMforError(uExInt,x,y,t,param);
    [u, ~, meshCoarse, ~, errTrue(iter), ~, errDet]  = Elliptic_Solver(2, el, vert, bound, 'P1', data);
    errTrue(iter) = errTrue(iter) / solNorm;
    N = size(meshCoarse.elements, 2);
    uInt = pdeInterpolant(meshCoarse.vertices, meshCoarse.elements, u);
    data.uexact = @(x,y,t,param) evaluateFEMforError(uInt,x,y,t,param);
    
    NVec(iter) = N;
    
    % Perturb vertices & compute probabilistic solutions
    idxInternal = setdiff(1:size(meshCoarse.vertices, 2), [meshCoarse.boundaries(1, :), meshCoarse.boundaries(2, :)]);
    NVertInt = length(idxInternal);
    hPerturbation = zeros(1, NVertInt);
    itPert = 1;
    for ii = idxInternal
        [~, Neighbours] = find(meshCoarse.elements(1:3, :) == ii);
        hPerturbation(itPert) = min(meshCoarse.h(Neighbours));
        itPert = itPert + 1;
    end
    
    errProbSqd = zeros(M, N);
    for jj = 1 : M
        MESH = meshCoarse;
        MESH.vertices(:, idxInternal) = MESH.vertices(:, idxInternal) + 5 * repmat(hPerturbation,2,1).^2 .* (rand(2, NVertInt) - 0.5);
        [U, ~, ~, ~, ~, ~, L2Loc] = Elliptic_Solver(2, MESH.elements, MESH.vertices, MESH.boundaries, 'P1', data, [], [], false);
        errProbSqd(jj, :) = L2Loc.^2 ./ meshCoarse.h';
    end
    errProbSqd = mean(errProbSqd);
    errProb = sqrt(errProbSqd);
           
    % Refine mesh
    errEst(iter) = sqrt(sum(errProbSqd)) / solNorm;
    display(['error estimate ', num2str(errEst(iter))])
    mark = [];
    
    if errEst(iter) < tol
        flag = false;
    else
        for i = 1 : N
            if errProb(i) > tol * solNorm / sqrt(N)
                mark = [mark i];
            else
            end
        end
        if length(mark) < 0.01 * N
            flag = false;
        end
        [vert, bound, el] = refinemesh(g, vert, bound, el, [], mark', 'regular');
    end
    
    if iter == 1
        fig = createFigure(W, H, 'enhanced',enhanced);
        FEMplot2D(meshCoarse);
        title({'initial mesh'}, 'interpreter', 'latex')
        xlabel('$x_1$', 'interpreter', 'latex')
        ylabel('$x_2$', 'interpreter', 'latex')
        axis square
        set(gca, 'innerpos',inPos);
        %export_fig(fig, '../../../Reports/DraftPDE_18/VERSION7/Figures/MeshAdapt2D_InitMesh.eps', '-nocrop', '-gray')

    elseif ~flag
        fig = createFigure(W, H, 'enhanced',enhanced);
        plotP0(meshCoarse, errDet / solNorm)
        zlim = caxis;
        title({'true error'}, 'interpreter', 'latex')
        colorbar off
        xlabel('$x_1$', 'interpreter', 'latex')
        ylabel('$x_2$', 'interpreter', 'latex')
        axis square
        set(gca, 'innerpos',inPos);
        %export_fig(fig, '../../../Reports/DraftPDE_18/VERSION7/Figures/MeshAdapt2D_TrueErr.eps', '-nocrop', '-painters')

        fig = createFigure(W, H, 'enhanced',enhanced);
        plotP0(meshCoarse, errProb' / solNorm)
        caxis(zlim)
        title({'error estimator'}, 'interpreter', 'latex')
        colorbar off
        xlabel('$x_1$', 'interpreter', 'latex')
        ylabel('$x_2$', 'interpreter', 'latex')
        axis square
        set(gca, 'innerpos',inPos);
        %export_fig(fig, '../../../Reports/DraftPDE_18/VERSION7/Figures/MeshAdapt2D_ErrEst.eps', '-nocrop', '-painters')
        
        fig = createFigure(W, H, 'enhanced',enhanced);
        FEMplot2D(meshCoarse);
        title({'final mesh'}, 'interpreter', 'latex')
        xlabel('$x_1$', 'interpreter', 'latex')
        ylabel('$x_2$', 'interpreter', 'latex')
        axis square
        set(gca, 'innerpos',inPos);
        %export_fig(fig, '../../../Reports/DraftPDE_18/VERSION7/Figures/MeshAdapt2D_FinalMesh.eps', '-nocrop', '-gray')
        
        fig = createFigure(W, H, 'enhanced',enhanced);
        loglog(NVec, errEst, 'ks-')
        hold on
        loglog(NVec, errTrue, 'k<-')
        loglog(NVec, tol * ones(size(NVec)), 'k--')
        legend({'estimate', 'true', 'tol'}, 'position', [0.1758 0.1512 0.3837 0.1865])
        title({'error convergence'}, 'interpreter', 'latex')
        xlabel('$N$', 'interpreter', 'latex')
        ylabel('error', 'interpreter', 'latex')
%         ylim([1e-3, 1])
        axis square
        set(gca, 'innerpos',inPos);
        %export_fig(fig, '../../../Reports/DraftPDE_18/VERSION7/Figures/MeshAdapt2D_ErrConv.eps', '-nocrop', '-gray')
    end
end
