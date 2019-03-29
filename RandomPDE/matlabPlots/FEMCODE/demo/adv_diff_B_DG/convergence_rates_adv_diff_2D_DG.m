%% CONVERGENCE TEST OF DG FOR ADVECTION-DIFFUSION PROBLEM OF TYPE
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

scal_diffusion = 0.0001;
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

%% CONVERGENCE TESTS FOR P1-DG ELEMENTS ON 'Nmax' ITERATIVELY REFINED MESHES
Nmax = 7;

options.dgnorm = true;

h1 = zeros(Nmax,1);
l2 = zeros(Nmax,1);
h1s = zeros(Nmax,1);
dg = zeros(Nmax,1);
dof = zeros(Nmax,1);

sol = cell(Nmax,1);
femspace = cell(Nmax,1);
mesh = cell(Nmax,1);
vps = cell(Nmax,1);
eqn = cell(Nmax,1);

for i=1:Nmax
        
    % construct mesh
    N = 2^(i-1);    
    thismesh = structured_mesh([0,1,0,1],N,struct('centre',true));
    thismesh.bdflag = 'dirichlet';
    
    % calculate solution
    [sol{i}, femspace{i}, mesh{i}, vps{i}, eqn{i}] =...
        elliptic_dg(thismesh, vp, options);
    
    % calculate error
    dof(i) = femspace{i}.ndof;
    [h1(i), l2(i), h1s(i), dg(i)] = ...
        get_DGerror_adv_diff_exact(mesh{i}, femspace{i}, vps{i}, eqn{i}, sol{i}, @ex_adv_diff_reac1_2D);
    
    fprintf('P1-FEM with #DOF = %d DONE.\n',dof(i));
end

%% TURN ABSOLUTE ERROR INTO RELATIVE ERROR
[ H1norm, L2norm, H1seminorm] = get_H1norm(mesh{i}, femspace{i}, sol{i});
h1 = h1 / H1norm;
l2 = l2 / L2norm;
h2s = h1s / H1seminorm;

%% PLOT L2-ERROR
figure;
loglog(dof(:,1),l2(:,1),'b-x', ...
    [10^2,10^3],(([10^2,10^3]).^(-1/2)).^2,'b--', ...
    [10^2,10^3],(([10^2,10^3]).^(-1/2)).^1.5,'b:');
title('Error in L2-norm');
xlabel('DOF');
ylabel('relative error');
legend('p1','O(h2)','O(h1.5)');

%% PLOT H1-ERROR
figure;
loglog(dof(:,1),h1(:,1),'b-x', ...
    [10^2,10^3],(([10^2,10^3]).^(-1/2)).^1,'b--');
title('Error in H1-norm');
xlabel('DOF');
ylabel('relative error');
legend('p1','O(h)');

%% PLOT DG-ERROR
figure;
loglog(dof(:,1),dg(:,1),'b-x', ...
    [10^2,10^3],(([10^2,10^3]).^(-1/2)).^1,'b--', ...
    [10^2,10^3],(([10^2,10^3]).^(-1/2)).^1.5,'b:');
title('Error in DG-norm');
xlabel('DOF');
ylabel('absolute error');
legend('p1','O(h)','O(h1.5)');