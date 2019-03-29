%% CONVERGENCE TEST FOR NONLINEAR MONOTONE PROBLEM
%
%   -div(A(grad u)) = f(x)
%
% on the domain Omega = [0,1]^2.

%% WARNING, THIS TEST TAKES ~5 MINUTES
clear vp sol;

%% DEFINE DATA
vp.A = @operator_mon2d; % operator
vp.dxiA = @operator_jac_mon2d; % Jacobi of operator
vp.f = @rhs_mon2d;

%% CONVERGENCE TESTS FOR P1-,P2-,P3-ELEMENTS ON 'Nmax' ITERATIVELY REFINED MESHES
Nmax = 6;

% tolerance and maximal number of iterations for Newton's method
options.tolerance = 1e-6;
options.maxIt = 25;

% possibility to set initial guess for Newton's method
% options.sol = zeros(NDOF,1);

options.verbose = true;
% vp.solver = 'agmg';

h1 = zeros(Nmax,3);
l2 = zeros(Nmax,3);
h1s = zeros(Nmax,3);
dof = zeros(Nmax,3);

for j=1:3
    vp.elemtype = ['p' num2str(j)];
    fprintf('Starting P%d-FEM simulations...\n',j);
for i=1:Nmax
    N = 2^(i-1);    
    thismesh = structured_mesh([0,1,0,1],N,struct('centre',true));
    thismesh.bdflag = 'dirichlet';
    [sol{i,j}, femspace{i,j}, mesh{i,j}, vps{i,j}] = poisson_mon(thismesh, vp, options);
    dof(i,j) = femspace{i,j}.ndof;
    [h1(i,j), l2(i,j), h1s(i,j)] = ...
        get_H1error_exact(mesh{i,j}, femspace{i,j}, sol{i,j}, @ex2_mon);
    fprintf('P%d-FEM with #DOF = %d DONE.\n',j,dof(i,j));
end
end

%% PLOT L2-ERROR
figure;
loglog(dof(:,1),l2(:,1),'b-x', ...
    dof(:,2),l2(:,2),'g-x', ...
    dof(:,3),l2(:,3),'r-x', ...
    [10^2,10^3],(([10^2,10^3]).^(-1/2)).^2,'b--', ...
    5*[10^2,10^3],((5*[10^2,10^3]).^(-1/2)).^3,'g--', ...
    10*[10^2,10^3],((10*[10^2,10^3]).^(-1/2)).^4,'r--');
title('Error in L2-norm');
xlabel('DOF');
ylabel('absolute error');
legend('p1','p2','p3','O(h2)','O(h3)', 'O(h4)');

%% PLOT H1-ERROR
figure;
loglog(dof(:,1),h1(:,1),'b-x', ...
    dof(:,2),h1(:,2),'g-x', ...
    dof(:,3),h1(:,3),'r-x', ...
    [10^2,10^3],(([10^2,10^3]).^(-1/2)).^1,'b--', ...
    5*[10^2,10^3],((5*[10^2,10^3]).^(-1/2)).^2,'g--', ...
    10*[10^2,10^3],((10*[10^2,10^3]).^(-1/2)).^3,'r--');
title('Error in H1-norm');
xlabel('DOF');
ylabel('absolute error');
legend('p1','p2','p3','O(h)','O(h2)', 'O(h3)');


