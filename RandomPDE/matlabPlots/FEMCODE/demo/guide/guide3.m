%% GUIDE 3 - Stokes problem
% Ondrej Budac, 2014
%
% Here we describe how to solve
% Stokes problem with different boundary conditions, finite elements, and
% other options, using the ANMC-FE package. 

%% Problem formulation
% $$ -\textnormal{div}\, T(u(x), p(x))
% = f(x) \qquad \textnormal{in\ }\Omega$$
%
% $$ T(u(x), p(x)) = \mu(\nabla u(x) + \nabla^T u(x)) -p(x)I $$
%
% $$ \textnormal{div}\, u(x) = 0 \qquad \textnormal{in\ }\Omega $$
%
% $$ u(x) = g_D(x) \qquad \textnormal{on\ } \Gamma_D, $$
%
% $$ T(u(x), p(x)) \cdot n = g_N(x) 
% \qquad \textnormal{on\ } \Gamma_N,$$
%
% where $u(x)$ is a velocity field, $p(x)$ is a pressure field, $T$ is the
% stress tensor, $\mu$ is a dynamic viscosity, $f(x)$ is a volumetric force
% field, $\Gamma_D$ is the Dirichlet portion of the boundary 
% $\partial\Omega$, $\Gamma_N$ is the Neumann portion of the boundary, 
% $g_D(x)$ is the Dirichlet trace and $g_N(x)$ is the Neumann trace.
%
% Since we only solve the incompressible problem at the moment, we can
% simplify it to its more common description
%
% $$ -\mu\Delta u(x) + \nabla p(x) = f(x) \qquad \textnormal{in\ }\Omega$$
%
% $$ \textnormal{div}\, u(x) = 0 \qquad \textnormal{in\ }\Omega $$
%
% $$ u(x) = g_D(x) \qquad \textnormal{on\ } \Gamma_D, $$
%
% $$ \mu \frac{\partial u(x)}{\partial n}  - p(x)n = g_N(x) 
% \qquad \textnormal{on\ } \Gamma_N,$$


%% Dirichlet B.C.
% We've already seen how simple meshes can be constructed. Let us define a
% structured mesh on a unit square:

mesh = structured_mesh([-1,1,-1,1], 20, struct('centre',true));

%%
% The last parameter ensures that any triangle has at most one boundary
% edge. It is one of the Stability conditions for Stokes systems.

%%
% We start with a Stokes problem on a square domain with zero
% Dirichlet boundary conditions, unit viscosity, and circular force field
%
% $$ -\Delta u(x) + \nabla p(x) = f(x) \qquad \textnormal{in\ }\Omega$$
%
% $$ \textnormal{div}\, u(x) = 0 \qquad \textnormal{in\ }\Omega $$
%
% $$ u(x) = 0 \qquad \textnormal{on\ } \partial\Omega, $$
%
% We assume that $f(x) = (-x_2, x_1)$. 
%
% We need to define two finite element spaces. Standard theory for Stokes
% says that there are stable pairs of finite elements as the Taylor-Hood
% elements P2-P1. We will consider only this pair in the following
% examples. Another stable implemented pairs are P3-P2 and P1b-P1.

fh = @(x)([-x(:,2), x(:,1)]);
vp = struct('a', 1, 'f', fh, ...
    'elemtype', 'p2', 'pelemtype', 'p1', 'bc', [0, 0]);
mesh.bdflag = 'dirichlet';

[usol, ufemspace, psol, pfemspace, mesh, vp] = stokes(mesh, vp);

%%
% Note that the pressure solution of pure Dirichlet problem is not unique
% and hence a representant with zero average is computed. 
%
% The stokes solver behaves similarly as the poisson solver, but it gives
% more output. Since we are finding velocity AND pressure, there are two
% discrete solutions (|usol|, |psol|) and their corresponding finite
% element spaces (|ufemspace|, |pfemspace|). Let us plot the solutions now

% Force field
bary = get_rc(mesh); % computes barycenters of elements
fbary = fh(bary);
quiver(bary(:,1), bary(:,2), fbary(:,1), fbary(:,2)); snapnow;

%%
%

% pressure solution
clf; simpplot_sol(mesh, psol); snapnow;

%%
%

% velocity solution as quiver plot
N = size(mesh.node,1); 
clf; quiver(mesh.node(:,1), mesh.node(:,2), usol(1:N, 1), usol(1:N, 2));
snapnow;

%%
%

% velocity components
clf; simpplot_sol(mesh, usol(:,1)); snapnow;
clf; simpplot_sol(mesh, usol(:,2)); snapnow;

%% 
% *Non-homogeneous Dirichlet B.C.*
% When presecribing non-homogeneous Dirichlet Boundary conditions, 
%
% $$ \Delta u(x) + \nabla p(x) = f(x) \qquad \textnormal{in\ }\Omega$$
%
% $$ \textnormal{div}\, u(x) = g_D(x) \qquad \textnormal{in\ }\Omega $$
%
% $$ u(x) = 0 \qquad \textnormal{on\ } \partial\Omega, $$
%

%one has to ensure the compatibility, that is:
%
% $$ \int_{x \in \partial\Omega} g_D(x) \cdot n\, \textnormal{d}s = 0.$$
%
% We now solve this problem with $f(x) \equiv (0,-1)$ and Pouiseille inflow
% on the left side of the domain and Pouiseille outflow on the bottom side 
% of the domain.

mesh = structured_mesh([0,1,0,1], 20, struct('centre',true)); 
gDh = @(x)([(1-x(:,1)).*x(:,2).*(1-x(:,2)), -(1-x(:,2)).*x(:,1).*(1-x(:,1))]);
vp = struct('a', 1, 'f', [0,-1], ...
    'elemtype', 'p2', 'pelemtype', 'p1', 'bc', {{gDh, [0,0]}});
mesh.bdflag = 'dirichlet';

[usol, ~, psol] = stokes(mesh, vp);

%%

% velocity solution as quiver plot
N = size(mesh.node,1); 
clf; quiver(mesh.node(:,1), mesh.node(:,2), usol(1:N, 1), usol(1:N, 2));
snapnow;

%%
%

% pressure solution
clf; simpplot_sol(mesh, psol); snapnow;


%% Neumann B.C.
% We come directly to non-homogeneous Neumann B.C. Let us solve the problem
%
% $$ -\mu\Delta u(x) + \nabla p(x) = f(x) \qquad \textnormal{in\ }\Omega$$
%
% $$ \textnormal{div}\, u(x) = 0 \qquad \textnormal{in\ }\Omega $$
%
% $$ \mu \frac{\partial u(x)}{\partial n}  - p(x)n = g_N(x) 
% \qquad \textnormal{on\ } \Gamma_N,$$
%
% with $\mu = 1$, $f(x) = (0, 0)$ and $g_N = (\sin(\pi(x_1+x_2)), 
% \cos(\pi(x_1-x_2)))$. 

mesh = structured_mesh([-1,1,-1,1], 20, struct('centre',true));
gNh = @(x)([sin(pi*(x(:,1)+x(:,2))), cos(pi*(x(:,1)-x(:,2))) ]);
vp = struct('a', 1, 'f', [0, 0], ...
    'elemtype', 'p2', 'pelemtype', 'p1', 'bc', {{[0,0] , gNh}});
mesh.bdflag = 'neumann';

[usol, ~, psol] = stokes(mesh, vp);

%%
% Note that the velocity solution of pure Dirichlet problem is not unique
% and hence a representant (with zero average of every component) is
% computed. The problem is solvable only if an appropriate compatiblity
% condition is true:
%
% $$\int_\Omega f_i(x) = \int_{\partial \Omega} g_{N_i}(x))$$
%
% for every $i = 1, \ldots, d$. Let us plot the solutions.


% velocity solution as quiver plot
N = size(mesh.node,1); 
clf; quiver(mesh.node(:,1), mesh.node(:,2), usol(1:N, 1), usol(1:N, 2));
snapnow;

%%
%

% pressure solution
clf; simpplot_sol(mesh, psol); snapnow;

%% Mixed Boundary conditions
% Let us solve the full mixed problem
%
% $$ -\mu\Delta u(x) + \nabla p(x) = f(x) \qquad \textnormal{in\ }\Omega$$
%
% $$ \textnormal{div}\, u(x) = 0 \qquad \textnormal{in\ }\Omega $$
%
% $$ u(x) = g_D(x) \qquad \textnormal{on\ } \Gamma_D, $$
%
% $$ \mu \frac{\partial u(x)}{\partial n}  - p(x)n = g_N(x) 
% \qquad \textnormal{on\ } \Gamma_N,$$
%
% where $\mu =1$, $f(x) \equiv (0, -10)$, $\Gamma_N$ is the right side of
% the square boundary and $\Gamma_D$ is the rest. We let $g_N(x) \equiv
% (0,1)$ and is the Pouseille inflow $g_D(x) = [(1-x_1)x_2(1-x_2), 0]$. We
% get the following problem description and solution.

gDh = @(x)([(1-x(:,1)).*x(:,2).*(1-x(:,2)), zeros(size(x,1),1)]);
vp = struct('a', 1, 'f', [0, -10], ...
    'elemtype', 'p2', 'pelemtype', 'p1', 'bc', {{gDh, [0,1]}});

mesh = structured_mesh([0,1,0,1], 20, struct('centre',true)); 
mesh.bdflag =cell(2,1);
mesh.bdflag{1} = get_bdflag(mesh, @(x)((x(:,1)<1-1e-5)));
mesh.bdflag{2} = get_bdflag(mesh, @(x)((x(:,1)>1-1e-5)));
[usol, ~, psol] = stokes(mesh, vp);

%%

% velocity solution as quiver plot
N = size(mesh.node,1); 
clf; quiver(mesh.node(:,1), mesh.node(:,2), usol(1:N, 1), usol(1:N, 2));
snapnow;

%%
%

% pressure solution
clf; simpplot_sol(mesh, psol); snapnow;