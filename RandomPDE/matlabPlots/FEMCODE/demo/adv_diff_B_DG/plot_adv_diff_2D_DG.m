%% PLOT DG-SOLUTION FOR ADVECTION-DIFFUSION PROBLEM OF TYPE
%
%   -div(a grad(u)) + b * grad(u) = f(x)
%
% on the domain Omega = [0,1]^2 with constant data a, b and exact
% solution u(x).

clear vp sol;

%% DEFINE THE DIFFUSION TENSOR, ADVECTION, RHS-FUNC, BOUNDARY CONDITION
% 1) scaling of the different terms (to produce diffusion or advection or
%    reaction dominated problems without changing the 'core' data
% 2) diffusion a given by scalar a1
% 3) advection b given by vector [b1;b2]
% 4) exact solution u(x) defined in ex_adv_diff_reac1.m
% 5) right-hand side f(x) automatically defined by 1-4
% 6) zero Dirichlet boundary conditions

scal_diffusion = 0.00001;
scal_advection = 1;

a1 = 1;
vp.a = scal_diffusion*a1;

b1 = 1; b2 = 0.5;
vp.b = scal_advection*[b1;b2];

vp.f = @(x,k) -(vp.a.*ex_adv_diff_reac1(x,[2,0])...
                + vp.a.*ex_adv_diff_reac1(x,[0,2]))...
              + vp.b(1).*ex_adv_diff_reac1(x,[1 0])...
              + vp.b(2).*ex_adv_diff_reac1(x,[0 1]);

%% P1-DG ELEMENTS (solving advection with DG only implemented for P1 elements)
vp.elemtype = 'p1d';

%% CHOOSE THE DG-FORM FOR THE DIFFUSION
% vp.diffpart = 'sipg'; % SYMMETRIC INTERIOR PENALTY
% vp.diffpart = 'nipg'; % NON-SYMMETRIC INTERIOR PENALTY
vp.diffpart = 'iipg'; % INCOMPLETE INTERIOR PENALTY

% penalization parameter for diffusive part
vp.alpha = 10;

%% CHOOSE SUFFICIENT QUADRATURE RULES FOR RHS
% diffusion dominated
if scal_diffusion >= scal_advection
    vp.fquad = 1;
else
    % advection dominated
    vp.fquad = 2;
end

options.verbose = true;

%% CONSTRUCT MESH AND CALCULATE SOLUTION
refine = 6;
N = 2^(refine-1);    
thismesh = structured_mesh([0,1,0,1],N,struct('centre',true));
thismesh.bdflag = 'dirichlet';
[sol, femspace, mesh, vps] = elliptic_dg(thismesh, vp, options);

%% PLOT DG-SOLUTION
simpplot_sol_dg(mesh,sol,femspace);