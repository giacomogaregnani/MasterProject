%% GUIDE 2 - Elliptic problem
% Ondrej Budac, 2014
%
% Here we describe how to solve
% Poisson problem with different boundary conditions, finite elements, and
% other options, using the ANMC-FE package. 

%% Problem formulation
% We consider only an elliptic problem with no advection nor
% convection. Fo these phenomena, see another guide. We have
%
% $$ -\textnormal{div}(a(x)\nabla u(x))
% = f(x) - \textnormal{div}(a(x)g(x))\qquad \textnormal{in\ }\Omega,$$
%
% $$ u(x) = g_D(x) \qquad \textnormal{on\ } \Gamma_D, $$
%
% $$ a(x)(\nabla u(x) - g(x)) \cdot n 
% = g_N(x) \qquad \textnormal{on\ } \Gamma_N, $$
%
% where $f(x)$ is a scalar source, $g(x)$ is a vector field, 
% $a(x)$ is a positive
% definite tensor field, boundary $\partial \Omega$ is divided into a
% Dirichlet part ($\Gamma_D$) and a Neumann part ($\Gamma_N$), and scalar
% functions $g_D$ and $g_N$ define the Dirichlet and Neumann boundary
% traces, respectively.
%
% We start with a simpler problem and then we introduce generalization
% leading to this form.

%% First problem
% We've already seen how simple meshes can be constructed. Let us define a
% structured mesh on a unit square:

mesh = structured_mesh([0,1,0,1], 30); % shortcut for [30, 30]

%%
% We start with a Poisson problem on a square domain with zero
% Dirichlet boundary conditions and constant source term:
%
% $$ - \Delta u(x) = 1\qquad \textnormal{in\ }\Omega $$
%
% $$ u(x)  = 0\qquad \textnormal{on\ }\partial\Omega $$
%
% Here, the tensor $a(x)$ is an identity matrix, $a(x) \equiv I$ and the
% source term $f(x) \equiv 1$ is constant.
%
% This problem can be defined as follows. 

vp = struct('a', eye(2), 'f', 1, 'elemtype', 'p1');
mesh.bdflag = 'dirichlet';

%%
% We specify Dirichlet boundary conditions in the mesh and they are 
% homogeneous (default). We use P1 finite elements and a constant 
% right-hand side. The Solution of this discrete problem is then found via

[sol, femspace, mesh, vp] = poisson(mesh, vp);

%%
% Before explaining what are the outputs, let us plot the solution using
% the following routine.

simpplot_sol(mesh,sol); snapnow; 

%%
% There are four outputs of the |poisson| script. Let us describe them 
% in the reversed order
%
% * |vp| - the solver will fill values of the variational problem (e.g.,
% quadrature formulas, some flags, etc.) and it will return this updated 
% structure.

vp

%%
% * |mesh| - the same goes for the mesh, help fields as dlambda, bdflag, or
% volume will appear or change.

mesh

%%
% * |femspace| - this is a description of the finite element spaces.
% Finite element spaces are also structures and they posses some
% information about how degrees of freedom are mapped into the mesh.
% This femspace contains its elemtype ('p1'), number of degrees of freedom
% (ndof), local number of degrees of freedom in one element (ldof). its
% polynomial degree (deg) and two other fields (bd_dof contains the
% boundary degrees of freedom, and nodelambda contains local barycentric
% coordinates of the nodes).

femspace

%%
% * |sol| - Finally, having a mesh and a femspace, one can encode a finite 
% element function in a single vector. This vector now has length of
% |femspace.ndof|

size(sol)

%% Equation
% Let us discover the scope of equations we can solve. 
%
% *Non-symmetric tensor*
%
% $$ - \textnormal{div}(a \nabla u(x)) = 1\qquad \textnormal{in\ }\Omega $$
%
% $$ u(x)  = 0\qquad \textnormal{on\ }\partial\Omega $$
%
% with $a = {2\ \ 1.8 \choose 0\quad 2}$.

vp = struct('a', [2,1.8; 0,2], 'f', 1, 'elemtype', 'p1');

mesh = structured_mesh([0,1,0,1], 30); 
mesh.bdflag = 'dirichlet';

sol = poisson(mesh, vp);

% plot
clf; simpplot_sol(mesh,sol); snapnow; 

%%
% *Varying tensor*
%
% Let us vary the tensor $a(x)$ to get the problem
%
% $$ -\textnormal{div}(a(x)\nabla u(x))=1\qquad \textnormal{in\ }\Omega $$
%
% $$ u(x) = 0\qquad \textnormal{on\ }\partial\Omega $$
%
% This can be obtained programmatically by passing a function handle in the
% field |vp.a|. If
%
% $$ a_{11} = 0.1+x_1^2, \qquad a_{22} = 0.1 + x_2^2, \qquad a_{12} = a_{21} = 0$$
%
% Then we can type

vp = struct('a', @guide2_a1, 'f', 1, 'elemtype', 'p1'); 
mesh.bdflag = 'dirichlet';
sol = poisson(mesh, vp);

dbtype guide2_a1.m
clf; simpplot_sol(mesh,sol); snapnow; 

%%
% We assume that the handle in |vp.a| takes three parameters: 1) an array
% of points $x$ where it should be computed, 2+3) coordinates $k$ and $l$ 
% of the tensor to return: $a_{kl}(x)$.

%%
% *Discrete tensor description (Advanced)*
%
% Since the integration of the weak formulation can't be done exactly, we
% use a quadrature formula. If computation of tensor values $a(x)$ 
% is expensive, one can
% preevaluate them on the quadrature points of FEM and pass it to the
% solver. Then vp.a is an array of size |NT x d x d x NQ|, where |NT|
% is the number of triangles in the mesh and |NQ| is the number of
% quadrature points inside every triangle. To reconstruct the previous
% example, we could call:

vp = struct('f', 1, 'elemtype', 'p1');
mesh.bdflag = 'dirichlet';

d = size(mesh.node,2);       % dimension = 2
[lambda, weight] = quadpts(d,1); % quadrature formula
NQ = numel(weight);                % number of quadrature points
NT = size(mesh.elem,1);            % number of triangles
vp.a = zeros(NT,d,d,NQ);       % initialize
for i=1:NQ
    % compute quadrature points
    x = get_rc(mesh, [], 'all', [], lambda(i,:));
    for k=1:d
        for l=1:d
            vp.a(:,k,l,i) = guide2_a1(x, k, l);
        end
    end
end

%%
% now we should get the same result as before:

sol2 = poisson(mesh, vp);
max(abs(sol-sol2)) % test the $L^\infty$ difference

%%
% *Right-hand side in H-1*
%
% If one wants to solve
%
% $$ - \textnormal{div}(a(x)\nabla u(x)) 
% = f - \textnormal{div}(a(x)g)\qquad \textnormal{in\ }\Omega$$
%
% for a given vector $g$ (say $g=(0,1)^T$) then we simply define the RHS vector
% |vp.f| with additional terms (columns |2| to |d+1| describing $g_1$ to 
% $g_d$):

                                % [f, g_1, g_2]
vp = struct('a', @guide2_a1, 'f', [1, 1,   1], 'elemtype', 'p1'); 
mesh.bdflag = 'dirichlet';
sol = poisson(mesh, vp);

clf; simpplot_sol(mesh,sol); snapnow; 

%%
% *Varying right-hand side*
%
% To describe varying RHS (whether in L2 or in H-1) is as simple as to 
% describe a varying tensor. Solving
%
% $$ - \textnormal{div}(a(x)\nabla u(x)) 
% = 1 - \textnormal{div}(a(x)g)\qquad \textnormal{in\ }\Omega$$
%
% with Dirichlet BC and
% with $f(x) = 2(1-x_1)^2$, $g_1(x) = -x_2$, and $g_2(x) = x_1$ is then
% performed by we simply type

vp = struct('a', eye(2),  'f', @guide2_f1, 'elemtype', 'p1'); 
mesh.bdflag = 'dirichlet';
sol = poisson(mesh, vp);

dbtype guide2_f1.m
clf; simpplot_sol(mesh,sol); snapnow; 

%%
% *Discrete right-hand side description (Advanced)*
%
% Since the integration of the weak formulation can't be done exactly, we
% use a quadrature formula. If computation of the values $f(x)$ and/or
% $g(x)$ is expensive, one can
% preevaluate them on the quadrature points of FEM and pass it to the
% solver. Then |vp.f| is an array of size |NT x (1 or (d+1)) x 1 x NQ|, 
% where |NT|
% is the number of triangles in the mesh and |NQ| is the number of
% quadrature points inside every triangle. To reconstruct the previous
% example, we could call:

vp = struct('a', eye(2), 'elemtype', 'p1');
mesh.bdflag = 'dirichlet';

d = size(mesh.node,2);       % dimension = 2
[lambda, weight] = quadpts(d,1); % quadrature formula
NQ = numel(weight);                % number of quadrature points
NT = size(mesh.elem,1);            % number of triangles
vp.f = zeros(NT,d+1,1,NQ);      % initialize
for i=1:NQ
    % compute quadrature points
    x = get_rc(mesh, [], 'all', [], lambda(i,:));
    rhs = guide2_f1(x);
    for k=1:d+1        
        vp.f(:,k,1,i) = rhs(:,k);
    end
end

%%
% now we should get the same result as before:

sol2 = poisson(mesh, vp);
max(abs(sol-sol2)) % test the $L^\infty$ difference


%% Changing the discretization
%
% To use e.g. P2 finite elements, we simply type

vp = struct('a', eye(2), 'f', 1, 'elemtype', 'p2');
mesh.bdflag = 'dirichlet';

sol = poisson(mesh, vp);
clf; simpplot_sol(mesh,sol); snapnow; 

%%
% The plotting routine still works since (as it is), it ignores the higher
% degree corrections and displays only a piece-wise linear plot.
%
% There are four continuous element types implemented: p1, p2, p3,
% and p1b, that is p1 elements with bubbles.
%
% There are discontinuous elements as p0, p1d, p2d, p3d, which we won't
% discuss here.

%% Boundary conditions
% Until now we only used the Dirichlet boundary conditions.
%
% *Periodic BC*
%
% To solve the problem
%
% $$ - \textnormal{div}(a(x)\nabla u(x)) = f(x) - \textnormal{div}(a(x)g(x))\qquad \textnormal{in\ }\Omega$$
%
% on a torus we can do the following. Let us, for simplicity, take $f = 0$,
% $g = (1, 0)$, and $a = (1+0.45(\sin(2\pi x_1)+\sin(2\pi x_2)))I$.
%
% We now create a periodic mesh (which imposes the periodic boundary
% conditions) but we still leave the field |bdflag| in mesh to be
% |dirichlet| (here it is not important, but there can be some further
% inner boundary in a problem).

mesh_per = structured_mesh([0,1,0,1], 15, struct('periodic', true)); 
mesh_per.bdflag = 'dirichlet';
vp = struct('a', @guide2_a2, 'f', [0, 1, 0], 'elemtype', 'p2'); 
sol = poisson(mesh_per, vp);

dbtype guide2_a2.m
clf; simpplot_sol(mesh_per,sol); snapnow; 

%%
% *Non-Homogeneous Dirichlet BC*
%
% To solve the problem
%
% $$ - \textnormal{div}(a(x)\nabla u(x)) = f(x) - \textnormal{div}(a(x)g(x))\qquad \textnormal{in\ }\Omega$$
%
% $$ u(x) = g_D(x) \qquad \textnormal{on\ } \partial\Omega.$$
%
% with, for example, $g_D = \sin(3x_1 + x_2)$, we can do the following. The
% structure |vp| needs another field |bc| that describes the Dirichlet
% (|vp.bc{1}|) and Neumann (|vp.bc{2}|) boundary conditions. Since we want
% only Dirichlet B.C., we set the second cell as |0|.

vp = struct('a', eye(2), 'f', 1, 'elemtype', 'p1', 'bc', {{@guide2_gD1, 0}});
mesh.bdflag = 'dirichlet';
sol = poisson(mesh, vp);

dbtype guide2_gD1.m
clf; simpplot_sol(mesh, sol); snapnow;

%%
% *NATURAL BC*
%
% To solve the problem 
%
% $$ - \textnormal{div}(a(x)\nabla u(x)) = f(x) - \textnormal{div}(a(x)g(x))\qquad \textnormal{in\ }\Omega$$
%
% $$ - (a(x)\nabla u(x) - g(x)) \cdot n = g_N(x) \qquad \textnormal{on\ } \partial\Omega. $$
%
% with for example $g_N(x) = \nabla (\sin(x)\cos(y))\cdot n$ we can simply
% type the following.

vp = struct('a', eye(2), 'f', 0, 'elemtype', 'p1', 'bc', {{0, @guide2_gN1}});
mesh.bdflag = 'neumann';
sol = poisson(mesh, vp);

dbtype guide2_gN1.m
clf; simpplot_sol(mesh, sol); snapnow;

%%
% *MIXED BC*
%
% Consider the problem:
%
% $$ - \textnormal{div}(a(x)\nabla u(x)) = f(x) - \textnormal{div}(a(x)g(x))\qquad \textnormal{in\ }\Omega$$
%
% $$ u(x) = g_D(x) \qquad \textnormal{on\ } \Gamma_D. $$
%
% $$ - (a(x)\nabla u(x) - g(x)) \cdot n = g_N(x) \qquad \textnormal{on\ } \Gamma_N. $$
%
% We first specify the different boundaries in the field |mesh.bdflag{1}|
% (Dirichlet) and |mesh.bdflag{2}| (Neumann). We can do that using the
% function |get_bdflag| that accepts a mesh and a function handle (giving
% true or false, according to whether a point IS or ISN'T in the desired
% part of the boundary).
%
% Let $\Gamma_D$ be the left and right side of the squared domain and let
% $\Gamma_N$ be the upper and lower side. Suppose that $a$ is a unit
% matrix, $f(x) \equiv 0$, $g_D(x) \equiv 1$, and $g_N(x) \equiv -1$. Such
% a problem can be described as follows.

vp = struct('a', eye(2), 'f', 0, 'elemtype', 'p1', 'bc', {{1, -1}});
mesh.bdflag =cell(2,1);
mesh.bdflag{1} = get_bdflag(mesh, @(x)((x(:,1)<1e-5) | (1-x(:,1)<1e-5)));
mesh.bdflag{2} = get_bdflag(mesh, @(x)((x(:,2)<1e-5) | (1-x(:,2)<1e-5)));

[sol, femspace] = poisson(mesh, vp);
clf; simpplot_sol(mesh, sol); snapnow;
