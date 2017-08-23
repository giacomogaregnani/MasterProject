clc; clear; close all;

%% INVERSE PROBLEM 1D TEST (CONTINUOUS KAPPA - K.L. Expansion)

xKL = linspace(0, 1, 1000);
gamma = 0.01;
nu = 1;

C = @(x, y) nu * exp(-(x - y).^2 / gamma);

[XX, ~] = meshgrid(xKL, xKL);
CMatrix = C(XX, XX');

nEig = 5;
[V, D] = eigs(CMatrix, nEig);
D = abs(D);

clear XX CMatrix

%% GENERATE OBSERVATIONS

% Space discretization
meshRef = struct('xMin', 0, 'xMax', 1);
meshRef.N = 10000;
meshRef = buildMesh(meshRef);

% True field
% thetaRef = random('normal', zeros(nEig, 1), ones(nEig, 1));
% kappaRef = buildField_Continuous(thetaRef, D, V, xKL);
kappaRef = @(x) 1 + 2 * (x <= 0.75) .* (x >= 0.25);

% Force
% f = @(x) sin(2 * pi * x)';
f = @(x) 4 * x';

% Right boundary condition
rBC = 0;

% Solve
uRef = solveFwdProblemProb_Cont(meshRef, kappaRef, f, rBC);

% Generate observations
nObs = 10;
xObs = linspace(0.05, 0.95, nObs);
obsNoise = 10^(-4);

uObs = interp1(meshRef.x, uRef, xObs)' + obsNoise * randn(nObs, 1);

%% MCMC

% Runtime mesh
mesh = struct('xMin', 0, 'xMax', 1);

% prior choice
prior = struct('avg', [], 'stddev', []);
prior.stddev = eye(nEig);
prior.avg = zeros(nEig, 1);
prior.distr = 'normal';

plotDensities = true;

for i = 40
    nMCMC = 1e4;
    mesh.N = i; mesh = buildMesh(mesh);
    
    prob = true;
    [thetaAllProb, accRatioProb, SProb] = InverseProblemPDE1D_MCMC_Continuous(mesh, prior, f, rBC, nMCMC, prob, xObs, uObs, obsNoise, ...
                                                                             D, V, xKL);
    thetaAllProbEff = thetaAllProb(:, ceil(nMCMC / 10) : end);
    [pTprob, pUprob] = plotMCMCresultsPDE_Continuous(meshRef, thetaAllProbEff, kappaRef, f, rBC, uRef, D, V, xKL, xObs, uObs);
    yLim = get(pTprob.CurrentAxes, 'yLim');
    yLimU = get(pUprob.CurrentAxes, 'yLim');
    
    prob = false;
    [thetaAll, accRatio, S] = InverseProblemPDE1D_MCMC_Continuous(mesh, prior, f, rBC, nMCMC, prob, xObs, uObs, obsNoise, ...
                                                                  D, V, xKL);
    thetaAllEff = thetaAll(:, ceil(nMCMC / 10) : end);
    [pTdet, pUdet] = plotMCMCresultsPDE_Continuous(meshRef, thetaAllEff, kappaTheta, f, rBC, uRef, D, V, xKL, xObs, uObs);
    
    for j = 1 : 4
        figure
        hold on
        plotDensitiesMCMC(thetaAllEff(j, :));%, thetaRef(j));
        plotDensitiesMCMC(thetaAllProbEff(j, :));%, thetaRef(j));
    end
    
%     close all      
%     save(['resultsKL_', num2str(i)])
end

%

%% Additional plots

contours = false;
thetaTest = thetaAll;

if contours
    figure
    n1 = 5; n2 = 1;
    thetaGrid1 = linspace(min(thetaTest(n1, :)), max(thetaTest(n1, :)), 500);
    thetaGrid2 = linspace(min(thetaTest(n2, :)), max(thetaTest(n2, :)), 500);
    dens = akde(thetaTest([n1,n2], :)', thetaTest([n1,n2], :)');
    [XX, YY] = meshgrid(thetaGrid1, thetaGrid2);
    ZZ = griddata(thetaTest(n1, :), thetaTest(n2, :), dens, XX, YY);
    ZZ(isnan(ZZ)) = 0;
    contourf(XX, YY, ZZ, 12)
    hold on
    plot(thetaRef(n1), thetaRef(n2), '.r', 'markersize', 40)
end

scatters = false;

if scatters
    for i = 1 : 8
        figure
        %         scatter(thetaAllProbEff(i, :), thetaAllProbEff(i+1, :), 0.5, '.')
        hold on
        scatter(thetaAllEff(i, :), thetaAllEff(i+1, :), 0.5, '.r')
        plot(thetaRef(i), thetaRef(i+1), '.k', 'markersize', 30)
    end
end



