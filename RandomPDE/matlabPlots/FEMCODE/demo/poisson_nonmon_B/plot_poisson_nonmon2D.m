%% PLOT SOLUTION FOR NONLINEAR NONMONOTONE PROBLEM
%
%   -div(a(u) grad u) = f(x)
%
% on the domain Omega = [0,1]^2.

% TESTCASE OF ARTICLE
% A. Abdulle, G. Vilmart
% A priori error estimates for finite element methods with numerical
% quadrature for nonmonotone nonlinear elliptic problems, to appear in
% Numerische Mathematik

clear vp sol;

%% DEFINE DATA
vp.a = @tensor_nonmon2d; % tensor
vp.dsa = @tensor_jac_nonmon2d; % derivative of tensor
vp.f = @rhs_nonmon2d;

%% CHOOSE FEM-SPACE AND SETUP FOR NEWTON ITERATIONS
vp.elemtype = 'p1';

% tolerance and maximal number of iterations for Newton's method
options.tolerance = 1e-6;
options.maxIt = 100;

% possibility to set initial guess for Newton's method
% options.sol = zeros(NDOF,1);

options.verbose = 'true';

%% BUILD MESH AND COMPUTE SOLUTION
refine = 3;
N = 2^(refine-1);    
thismesh = structured_mesh([0,1,0,1],N,struct('centre',true));
thismesh.bdflag = 'dirichlet';
[sol, femspace, mesh, vps] = poisson_nonmon(thismesh, vp, options);

%% PLOT SOLUTION
simpplot_sol(mesh,sol,femspace);