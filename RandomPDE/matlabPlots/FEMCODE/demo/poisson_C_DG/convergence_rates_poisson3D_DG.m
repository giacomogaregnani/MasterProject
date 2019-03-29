%% CONVERGENCE TEST OF DG FOR DIFFUSION PROBLEM OF TYPE
%
%   -div(a grad(u)) = f(x)
%
% on the domain Omega = [0,1]^3 with constant data a, b and exact
% solution u(x).

%% WARNING: THIS TEST TAKES ~5 MINUTES
clear vp sol;

%% DEFINE DATA
vp.a = 1;
vp.f = @(x,k)(3*pi^2* sin(pi*x(:,1)) .* sin(pi*x(:,2)) .* sin(pi*x(:,3)) );

%% CHOOSE THE DG-FORM FOR THE DIFFUSION
vp.diffpart = 'sipg'; % SYMMETRIC INTERIOR PENALTY
% vp.diffpart = 'nipg'; % NON-SYMMETRIC INTERIOR PENALTY
% vp.diffpart = 'iipg'; % INCOMPLETE INTERIOR PENALTY

% penalization parameter
vp.alpha = 10;

%% CONVERGENCE TESTS FOR P1-,P2-DG ELEMENTS ON 'Nmax' ITERATIVELY REFINED MESHES
Nmax = 5;

vp.solver='agmg';

h1 = zeros(Nmax,2);
l2 = zeros(Nmax,2);
h1s = zeros(Nmax,2);
dg = zeros(Nmax,2);
dgs = zeros(Nmax,2);
dof = zeros(Nmax,2);

options.verbose = true;

for j=1:2
    vp.elemtype = ['p' num2str(j) 'd'];
    fprintf('Starting P%d-FEM simulations...\n',j);
for i=1:Nmax
    N = 2^(i-1);    
    thismesh = structured_mesh([0,1,0,1,0,1],N,struct('centre',true));
    thismesh.bdflag = 'dirichlet';
    [sol{i,j}, femspace{i,j}, mesh{i,j}, vps{i,j}] = poisson_dg(thismesh, vp, options);
    dof(i,j) = femspace{i,j}.ndof;
    [h1(i,j), l2(i,j), h1s(i,j), dg(i,j), dgs(i,j)] = ...
        get_DGerror_diffusion_exact(mesh{i,j}, femspace{i,j}, vps{i,j}, sol{i,j}, @ex3_dg);
    fprintf('P%d-FEM with #DOF = %d DONE.\n',j,dof(i,j));
end
end

%% TURN ABSOLUTE ERROR INTO RELATIVE ERROR
[ H1norm, L2norm, H1seminorm] = get_H1norm(mesh{i}, femspace{i}, sol{i});
h1 = h1 / H1norm;
l2 = l2 / L2norm;
h2s = h1s / H1seminorm;

%% PLOT L2-ERROR
figure;
loglog(dof(:,1),l2(:,1),'b-x', ...
    dof(:,2),l2(:,2),'g-x', ...
    [10^3,10^4],((0.2*[10^2,10^3]).^(-1/3)).^2,'b--', ...
    5*[10^3,10^4],((5*[10^2,10^3]).^(-1/3)).^3,'g--');
title('Error in L2-norm');
xlabel('DOF');
ylabel('relative error');
legend('p1','p2','O(h2)','O(h3)');

%% PLOT H1-ERROR
figure;
loglog(dof(:,1),h1(:,1),'b-x', ...
    dof(:,2),h1(:,2),'g-x', ...
    [10^3,10^4],((0.2*[10^2,10^3]).^(-1/3)).^1,'b--', ...
    5*[10^3,10^4],(([10^2,10^3]).^(-1/3)).^2,'g--');
title('Error in H1-norm');
xlabel('DOF');
ylabel('relative error');
legend('p1','p2','O(h)','O(h2)');

%% PLOT DG-ERROR
figure;
loglog(dof(:,1),dg(:,1),'b-x', ...
    dof(:,2),dg(:,2),'g-x', ...
    [10^3,10^4],((0.2*[10^2,10^3]).^(-1/3)).^1,'b--', ...
    5*[10^3,10^4],(([10^2,10^3]).^(-1/3)).^2,'g--');
title('Error in DG-norm');
xlabel('DOF');
ylabel('absolute error');
legend('p1','p2','O(h)','O(h2)');