%% DEPENDENCE ON PECLET NUMBER OF DG FOR ADVECTION-DIFFUSION PROBLEM WITH BOUNDARY LAYER
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

%% FIXED MESH, VARYING PECLET NUMBER
steps = 10;  

l2 = zeros(steps+1,1);
pe = zeros(steps+1,1);

sol = cell(steps+1);
femspace = cell(steps+1);
mesh = cell(steps+1);
vps = cell(steps+1);

% build mesh
refine = 6;
N = 2^(refine-1);    
thismesh = structured_mesh([0,1,0,1],N,struct('centre',true));  
thismesh.bdflag = 'dirichlet';

for i=0:steps
    
    %% DEFINE THE DIFFUSION TENSOR, ADVECTION, RHS-FUNC, BOUNDARY CONDITION
    % 1) scaling of the different terms (to produce diffusion or advection or
    %    reaction dominated problems without changing the 'core' data
    % 2) diffusion a given by scalar a1
    % 3) advection b given by vector [b1;b2]
    % 4) exact solution u(x) defined in ex_adv_diff_reac1.m
    % 5) right-hand side f(x) automatically defined by 1-4
    % 6) zero Dirichlet boundary conditions
    scal_diffusion = 1*10^(-i);
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
    
    %% CHOOSE SUFFICIENT QUADRATURE RULES FOR RHS
    % diffusion dominated
    if scal_diffusion >= scal_advection
        vp.fquad = 1;
    else
        % advection dominated
        vp.fquad = 2;
    end

    %% CALCULATE SOLUTION
    [sol{i+1}, femspace{i+1}, mesh{i+1}, vps{i+1}] = elliptic_dg(thismesh, vp);
    
    %% CALCULATE L2 ERROR
    pe(i+1) = scal_advection/scal_diffusion;
    [l2(i+1)] = ...
        get_L2error_exact(mesh{i+1}, femspace{i+1}, sol{i+1}, @(x,d)(ex_adv_diff_boundary_layer(x,scal_diffusion)));
    
    fprintf('P1-FEM with Pe = 1e%d DONE.\n',i);
    
end

%% PLOT L2-ERROR
figure;
loglog(pe,l2(:,1),'b-x', ...
    [pe(1),pe(3)],10e-4*[pe(1),pe(3)],'b--');
title('Error in L2-norm');
xlabel('pe');
ylabel('absolute error');
legend('p1','O(pe)');
