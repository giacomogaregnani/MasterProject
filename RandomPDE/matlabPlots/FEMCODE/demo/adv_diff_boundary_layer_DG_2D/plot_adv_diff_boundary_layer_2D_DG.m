%% GRAPHICAL COMPARISON DG- AND STANDARD-FEM-SOLUTION FOR BOUNDARY LAYER PROBLEM
%
%   -div(a grad(u)) + b * grad(u) = f(x)
%
% on the domain Omega = [0,1]^2 with constant data a, b and exact
% solution u(x).

clear vp sol;

%% P1-DG ELEMENTS (solving advection with DG only implemented for P1 elements)
vp.elemtype = 'p1d';

%% CHOOSE THE DG-FORM FOR THE DIFFUSION
% vp.diffpart = 'sipg'; % SYMMETRIC INTERIOR PENALTY
% vp.diffpart = 'nipg'; % NON-SYMMETRIC INTERIOR PENALTY
vp.diffpart = 'iipg'; % INCOMPLETE INTERIOR PENALTY

% penalization parameter for diffusive part
vp.alpha = 10;

%% DEFINE THE DIFFUSION TENSOR, ADVECTION, RHS-FUNC, BOUNDARY CONDITION
% 1) scaling of the different terms (to produce diffusion or advection or
%    reaction dominated problems without changing the 'core' data
% 2) diffusion a given by scalar a1
% 3) advection b given by vector [b1;b2]
% 4) exact solution u(x) defined in ex_adv_diff_reac1.m
% 5) right-hand side f(x) automatically defined by 1-4
% 6) zero Dirichlet boundary conditions
scal_diffusion = 0.01;
scal_advection = 1;

a1 = 1;
vp.a = scal_diffusion*a1;

b1 = 1; b2 = 1;
vp.b = scal_advection*[b1;b2];

a = 1-exp(-1/scal_diffusion);
b = @(x,y) exp(-(1-x).*(1-y)./scal_diffusion);
f_add = @(x,y) (x+y).*(1-b(x,y))/a + (2*y.*(1-y) + 2*x.*(1-x)).*b(x,y)./a + ...
    1/scal_diffusion*(x.*y.*(1-y).^2 + x.*y.*(1-x).^2 - x.*y.*(1-x) - x.*y.*(1-y)).*b(x,y)./a;
vp.f = @(x,k) f_add(x(:,1),x(:,2));

vp.bc = 'zero_dirichlet';

%% CHOOSE SUFFICIENT QUADRATURE RULES FOR RHS
% diffusion dominated
if scal_diffusion >= scal_advection
    vp.fquad = 1;
else
    % advection dominated
    vp.fquad = 2;
end

options.verbose = true;


%% BUILD MESH AND CALCULATE DG-SOLUTION
refine = 5;
N = 2^(refine-1);    
thismesh = structured_mesh([0,1,0,1],N,struct('centre',true));
[sol_dg, femspace_dg, mesh, vps_dg] = elliptic_dg(thismesh, vp, options);

%% PLOT DG-SOLUTION
simpplot_sol_dg(mesh,sol_dg,femspace_dg);

%% SETUP AND CALCULATE STANDARD-FEM SOLUTION
vp.elemtype = 'p1';

clear vp.alpha
clear vp.diffpart

[sol, femspace, mesh, vps] = elliptic(thismesh, vp, options);

%% PLOT STANDARD-FEM SOLUTION
options.view = 3;
options.edgecolor = 'black';
figure();
simpplot_sol(mesh,sol,options)