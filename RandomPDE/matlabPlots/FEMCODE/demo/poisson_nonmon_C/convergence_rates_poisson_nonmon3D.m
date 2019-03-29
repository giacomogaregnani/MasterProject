%% CONVERGENCE TEST FOR NONLINEAR NONMONOTONE PROBLEM
%
%   -div(a(u) grad u) = f(x)
%
% on the domain Omega = [0,1]^3.

% 3D EXTENSION OF TESTCASE OF ARTICLE
% A. Abdulle, G. Vilmart
% A priori error estimates for finite element methods with numerical
% quadrature for nonmonotone nonlinear elliptic problems, to appear in
% Numerische Mathematik

%% WARNING, THIS TEST TAKES ~19 MINUTES!
clear vp sol;

%% DEFINE DATA
vp.a = @tensor_nonmon3d; % tensor
vp.dsa = @tensor_jac_nonmon3d; % derivatives of tensor
vp.f = @rhs_nonmon3d;

%% CONVERGENCE TESTS FOR P1-,P2-ELEMENTS ON 'Nmax' ITERATIVELY REFINED MESHES
Nmax = 4;

% tolerance and maximal number of iterations for Newton's method
options.tolerance = 1e-6;
options.maxIt = 25;

% possibility to set initial guess for Newton's method
% options.sol = zeros(NDOF,1);

options.verbose = true;
vp.solver = 'agmg';

h1 = zeros(Nmax,2);
l2 = zeros(Nmax,2);
h1s = zeros(Nmax,2);
dof = zeros(Nmax,2);

for j=1:2
    vp.elemtype = ['p' num2str(j)];
    fprintf('Starting P%d-FEM simulations...\n',j);
for i=1:Nmax
    N = 2^(i-1);    
    thismesh = structured_mesh([0,1,0,1,0,1],N,struct('centre',true));
    thismesh.bdflag = 'dirichlet';
    [sol{i,j}, femspace{i,j}, mesh{i,j}, vps{i,j}] = poisson_nonmon(thismesh, vp, options);
    dof(i,j) = femspace{i,j}.ndof;
    [h1(i,j), l2(i,j), h1s(i,j)] = ...
        get_H1error_exact(mesh{i,j}, femspace{i,j}, sol{i,j}, @ex3_nonmon);
    fprintf('P%d-FEM with #DOF = %d DONE.\n',j,dof(i,j));
end
end

%% PLOT L2-ERROR
figure;
loglog(dof(:,1),l2(:,1),'b-x', ...
    dof(:,2),l2(:,2),'g-x', ...
    [10^2,10^3],(([10^2,10^3]).^(-1/3)).^2,'b--', ...
    5*[10^2,10^3],((5*[10^2,10^3]).^(-1/3)).^3,'g--');
title('Error in L2-norm');
xlabel('DOF');
ylabel('absolute error');
legend('p1','p2','O(h2)','O(h3)');

%% PLOT H1-ERROR
figure;
loglog(dof(:,1),h1(:,1),'b-x', ...
    dof(:,2),h1(:,2),'g-x', ...
    [10^2,10^3],(([10^2,10^3]).^(-1/3)).^1,'b--', ...
    5*[10^2,10^3],((5*[10^2,10^3]).^(-1/3)).^2,'g--');
title('Error in H1-norm');
xlabel('DOF');
ylabel('absolute error');
legend('p1','p2','O(h)','O(h2)');


