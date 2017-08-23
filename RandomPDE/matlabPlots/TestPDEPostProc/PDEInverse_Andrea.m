clc; clear; close all;

%% INVERSE PROBLEM 1D TEST

%% GENERATE OBSERVATIONS

% Space discretization
meshRef = struct('xMin', 0, 'xMax', 1);
meshRef.N = 10000;
meshRef = buildMesh(meshRef);

% True field
nParam = 2;
% thetaRef = random('normal', zeros(nParam, 1), ones(nParam, 1));
thetaRef = [0.5, -0.5, 0.5, -0.5]';
kappaRef = buildField_Andrea(thetaRef);

% Force
f = @(x) sin(2 * pi * x)';

% Right boundary condition
rBC = 0;

% data
data.f = f;
data.rBC = rBC;

% Solve
uRef = solveFwdProblem_Cont(meshRef, kappaRef, f, rBC);

% Generate observations
nObs = 3;
xObs = linspace(0.25, 0.75, nObs);
obsNoise = 10^(-3);

uObs = interp1(meshRef.x, uRef, xObs)' + obsNoise * randn(nObs, 1);

%% MCMC

% Runtime mesh
mesh = struct('xMin', 0, 'xMax', 1);

% prior choice
sPrior = eye(length(thetaRef));
mPrior = zeros(size(thetaRef));
prior = @(x) -0.5 * (x - mPrior)' * ((x - mPrior) ./ diag(sPrior));

% Initial guess
t0 = mPrior;

plotDensities = true;

for i = 5
    
    % Data
    nMCMC = 1e5;
    mesh.N = i; mesh = buildMesh(mesh);
    
    % Modeling error for corrected model
    nModError = 40;
    meshRefModErr = meshRef;
    [mE, vE] = estimateModErr(meshRefModErr, mesh, data, xObs, mPrior, sPrior, nModError);
    
    % Likelihood model
    likelihood = @(theta) likelihood_FEM(theta, xObs, uObs, obsNoise, mesh, data);
    likelihoodCorr = @(theta) likelihood_CorrFEM(theta, xObs, uObs, obsNoise, mE, vE, mesh, data);
    nMC = 10;
    likelihoodProb = @(theta) likelihood_ProbFEM(theta, xObs, uObs, obsNoise, mE, mesh, data, nMC);
    
    % Posteriors
    posterior = @(theta) prior(theta) + likelihood(theta);
    posteriorProb = @(theta) prior(theta) + likelihoodProb(theta);
    posteriorCorr = @(theta) prior(theta) + likelihoodCorr(theta);
    
    sigma = mesh.h;
    [thetaAll, accRatio, S] = InverseProblemPDE1D_MCMC_Andrea(t0, posterior, nMCMC, true, sigma);
    thetaAllEff = thetaAll(:, ceil(nMCMC / 10) : end);
    plotMCMCresultsPDE_Andrea(meshRef, thetaAllEff, thetaRef, f, rBC, uRef, xObs, uObs);
    
    [thetaAllCorr, accRatioCorr, SCorr] = InverseProblemPDE1D_MCMC_Andrea(t0, posteriorCorr, nMCMC, true, sigma);
    thetaAllCorrEff = thetaAllCorr(:, ceil(nMCMC / 10) : end);
    plotMCMCresultsPDE_Andrea(meshRef, thetaAllCorrEff, thetaRef, f, rBC, uRef, xObs, uObs);
    
    noisy = true;
    [thetaAllProb, accRatioProb] = InverseProblemPDE1D_MCMC_Andrea(t0, posteriorProb, nMCMC, true, sigma, noisy);
    thetaAllProbEff = thetaAllProb(:, ceil(nMCMC / 10) : end);
    plotMCMCresultsPDE_Andrea(meshRef, thetaAllProbEff, thetaRef, f, rBC, uRef, xObs, uObs);
    
end

%% Iterative prior

% prior choice
sPriorOut = eye(length(thetaRef));
mPriorOut = zeros(size(thetaRef)) - 1;
priorOut = @(x) -0.5 * (x - mPriorOut)' * ((x - mPriorOut) ./ diag(sPriorOut));
t0 = mPriorOut;

sPrior = sPriorOut;
mPrior = mPriorOut;

mesh = struct('xMin', 0, 'xMax', 1);
thetaGlobal = [];

nExp = 10;

for i = 1 : nExp
    
    % Data
    nMCMC = 5e3;
    mesh.N = 5; mesh = buildMesh(mesh);
    
    % Modeling error for corrected model
    nModError = 10;
    meshRefModErr = meshRef;
    [mE, vE] = estimateModErr(meshRefModErr, mesh, data, xObs, mPrior, sPrior, nModError);
    
    likelihoodCorr = @(theta) likelihood_CorrFEM(theta, xObs, uObs, obsNoise, mE, vE, mesh, data);
    posteriorCorr = @(theta) priorOut(theta) + likelihoodCorr(theta);
    
    sigma = mesh.h;
    if i == 1
        thetaNoChange = InverseProblemPDE1D_MCMC_Andrea(t0, posteriorCorr, nMCMC * nExp, true, sigma);
        plotMCMCresultsPDE_Andrea(meshRef, thetaNoChange, thetaRef, f, rBC, uRef, xObs, uObs);
    end
    
    [thetaAllCorr, accRatioCorr, SCorr] = InverseProblemPDE1D_MCMC_Andrea(t0, posteriorCorr, nMCMC, true, sigma);
    thetaGlobal = [thetaGlobal, thetaAllCorr];
    
    mPrior = mean(thetaAllCorr, 2);
    sPrior = cov(thetaAllCorr');
    prior = @(x) -0.5 * (x - mPrior)' * (sPrior \ (x - mPrior));
    sigma = SCorr;
    
    t0 = thetaAllCorr(:, end);
end

[pTdetCorr, pUdetCorr] = plotMCMCresultsPDE_Andrea(meshRef, thetaGlobal, thetaRef, f, rBC, uRef, xObs, uObs);



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