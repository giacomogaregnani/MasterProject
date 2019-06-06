%% PLOT FEM-SOLUTION FOR NONLINEAR MONOTONE PROBLEM
%
%   -div(A(grad u)) = f(x)
%
% on the domain Omega = [0,1]^2.

clear vp sol;

%% DEFINE DATA
vp.A = @operator_mon2d;
vp.dxiA = @operator_jac_mon2d;
vp.f = @rhs_mon2d;

%% CHOOSE FEM-SPACE AND SETUP FOR NEWTON ITERATIONS
vp.elemtype = 'p1';

% tolerance and maximal number of iterations for Newton's method
options.tolerance = 1e-6;
options.maxIt = 100;

% possibility to set initial guess for Newton's method
% options.sol = zeros(NDOF,1);

options.verbose = 'true';

%% BUILD MESH AND COMPUTE SOLUTION
refine = 4;
N = 2^(refine-1);    
thismesh = structured_mesh([0,1,0,1],N,struct('centre',true));
thismesh.bdflag = 'dirichlet';
[sol, femspace, mesh, vps] = poisson_mon(thismesh, vp, options);

%% PLOT SOLUTION
simpplot_sol(mesh,sol,femspace);