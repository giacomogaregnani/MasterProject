%% CONVERGENCE TEST FOR ADVECTION-DIFFUSION-REACTION PROBLEM OF TYPE
%
%   -div(a(x) grad(u)) + b(x) * grad(u) + c(x) u(x) = f(x)
%
% on the domain Omega = [0,1]^2 with general non-constant data a(x), b(x),
% c(x) and exact solution u(x) 

clear vp sol;

%% DEFINE THE DIFFUSION TENSOR, ADVECTION, REACTION, RHS-FUNC, BOUNDARY CONDITION
% 1) scaling of the different terms (to produce diffusion or advection or
%    reaction dominated problems without changing the 'core' data
% 2) diffusion a(x) defined in diffusion1.m
% 3) advection b(x) defined in advection1.m
% 4) reaction c(x) defined in reaction1.m
% 5) exact solution u(x) defined in ex_adv_diff_reac1.m
% 6) right-hand side f(x) automatically defined by 1-5
% 7) zero Dirichlet boundary conditions

scal_diffusion = 1;
scal_advection = 1;
scal_reaction = 1;

vp.a = @(x,k,l) scal_diffusion*diffusion2(x,k,l);

vp.b = @(x,k) scal_advection*advection2(x,k);

vp.c = @(x) scal_reaction*reaction2(x);


vp.D1a = @(x,k,l) scal_diffusion*diffusion2(x,k,l,[1 0]);
vp.D2a = @(x,k,l) scal_diffusion*diffusion2(x,k,l,[0 1]);

vp.f = @(x,k) -((vp.D1a(x,1,1) + vp.D2a(x,2,1)).*ex_adv_diff_reac2(x,[1 0])...
                + (vp.D1a(x,1,2) + vp.D2a(x,2,2)).*ex_adv_diff_reac2(x,[0,1])...
                + vp.a(x,1,1).*ex_adv_diff_reac2(x,[2,0])...
                + vp.a(x,2,2).*ex_adv_diff_reac2(x,[0,2])...
                + (vp.a(x,1,2) + vp.a(x,2,1)).*ex_adv_diff_reac2(x,[1,1]))...
              + vp.b(x,1).*ex_adv_diff_reac2(x,[1 0])...
              + vp.b(x,2).*ex_adv_diff_reac2(x,[0 1])...
              + vp.c(x).*ex_adv_diff_reac2(x);

options.verbose = false;

%% DEFINE ORDER OF QUADRATURE FORMULA APPLIED FOR DIFFUSION, ADVECTION AND REACTION
%
% OPTIONAL. Order can be different for the three terms

% vp.aquad = max([1,2*(femspace.deg - 1)]);
% vp.bquad = max([1,2*femspace.deg - 1]);
% vp.cquad = max([1,2*femspace.deg]);
% vp.fquad = max([1,2*(femspace.deg - 1)]);

%% CONVERGENCE TESTS FOR P1, P2, P3 ELEMENTS ON 'Nmax' ITERATIVELY REFINED MESHES
Nmax = 5;

for j=1:3
    vp.elemtype = ['p' num2str(j)];
    
    fprintf('Convergence test for %s elements started...\n',vp.elemtype);
    
    % DEFINE ORDER OF QUADRATURE FORMULA APPLIED FOR DIFFUSION, ADVECTION AND REACTION
    %
    % OPTIONAL. Order can be different for the three terms

    % DEFAULT CHOICE (INTERNAL CHOSEN IN THE ASSEMBLING)
    % vp.aquad = max([1,2*(j - 1)]);
    % vp.bquad = max([1,2*j - 1]);
    % vp.cquad = max([1,2*j]);
    % vp.fquad = max([1,2*(j - 1)]);
    
    % HOWEVER IT IS IMPORTANT TO AUGMENT THE ORDER OF THE QUADRATURE
    % FORMULA FOR THE RHS IF THE PROBLEM IS NOT DIFFUSION DOMINATED
   
    % diffusion dominated (optimal results, proofs available)
    vp.fquad = max([1,2*(j - 1)]);
    
    % advection dominated (not yet optimal results, influence of quadrature
    % error only studied for 'P1')
    % vp.fquad = max([1,2*j-1]);
    
    % reaction dominated (optimal results, but no proofs)
    % vp.fquad = max([1,2*j]);
    
for i=1:Nmax
    N = 2^(i-1);    
    thismesh = structured_mesh([0,1,0,1],N,struct('centre',true));
    thismesh.bdflag = 'dirichlet';
    [sol{i,j}, femspace{i,j}, mesh{i,j}, vps{i,j}] = elliptic(thismesh, vp, options);
    dof(i,j) = femspace{i,j}.ndof;
    [h1(i,j), l2(i,j), h1s(i,j)] = ...
        get_H1error_exact(mesh{i,j}, femspace{i,j}, sol{i,j}, @ex_adv_diff_reac2);
end
end

%% TURN ABSOLUTE ERROR INTO RELATIVE ERROR
[ H1norm, L2norm, H1seminorm] = get_H1norm(mesh{i,j}, femspace{i,j}, sol{i,j});
h1 = h1 / H1norm;
l2 = l2 / L2norm;
h2s = h1s / H1seminorm;

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
ylabel('relative error');
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
ylabel('relative error');
legend('p1','p2','p3','O(h)','O(h2)', 'O(h3)');
