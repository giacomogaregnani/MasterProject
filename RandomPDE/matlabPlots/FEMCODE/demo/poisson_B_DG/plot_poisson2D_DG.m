%% PLOT DG-SOLUTION FOR DIFFUSION PROBLEM OF TYPE
%
%   -div(a grad(u)) = f(x)
%
% on the domain Omega = [0,1]^2 with constant data a, b and exact
% solution u(x).

clear vp sol;

%% DEFINE DATA
vp.a = 1;
vp.f = @(x,k)(2*pi^2* sin(pi*x(:,1)) .* sin(pi*x(:,2)) );


%% CHOOSE THE DG-FORM FOR THE DIFFUSION
% vp.diffpart = 'sipg'; % SYMMETRIC INTERIOR PENALTY
% vp.diffpart = 'nipg'; % NON-SYMMETRIC INTERIOR PENALTY
vp.diffpart = 'iipg'; % INCOMPLETE INTERIOR PENALTY

% penalization parameter
vp.alpha = 10;

%% CHOOSE DG-SPACE (DG-PLOT ONLY FOR P1-DG)
vp.elemtype = 'p1d';

%% BUILD MESH AND CALC SOLUTION
refine = 5;
N = 2^(refine-1);    
thismesh = structured_mesh([0,1,0,1],N,struct('centre',true));
thismesh.bdflag = 'dirichlet';
[sol, femspace, mesh, vps] = poisson_dg(thismesh, vp);

%% PLOT DG SOLUTION
simpplot_sol_dg(mesh,sol,femspace);