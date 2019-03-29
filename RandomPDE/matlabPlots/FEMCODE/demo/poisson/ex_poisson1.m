function ex_poisson1
Nmax = 8;

h1 = zeros(Nmax,3);
l2 = zeros(Nmax,3);
h1s = zeros(Nmax,3);
dof = zeros(Nmax,3);

for j=1:3
	mesh = structured_mesh([0,1], 2);

	clear vp;
	vp.a = 1;
	vp.f = @rhs;
	vp.bc = 'zero_dirichlet';
	vp.elemtype = ['p' num2str(j)];
	for i=1:Nmax
		mesh = uniformrefine(mesh);
		[sol, femspace, mesh, vp] = poisson(mesh, vp);
		dof(i,j) = femspace.ndof;
		[h1(i,j), l2(i,j), ~] = ...
			get_H1error_exact(mesh, femspace, sol, @exactu);
	end
end

%%
[ H1norm, L2norm, H1seminorm] = get_H1norm(mesh, femspace, sol);
h1 = h1 / H1norm;
l2 = l2 / L2norm;
h2s = h1s / H1seminorm;

%%
figure;
loglog(dof(:,1),h1(:,1),'-x', ...
    dof(:,2),h1(:,2),'-x', ...
    dof(:,3),h1(:,3),'-x', ...
    [10^1,10^2],[10^-1,10^-2], ...
    5*[10^1,10^2],[10^-3,10^-5], ...
    10*[10^1,10^2],[10^-5,10^-8]/3);
legend('p1','p2','p3','O(h)','O(h2)', 'O(h3)');
end

function xx = exactu(x,der) %exactu exact solution
if nargin<2 || (der==0)
    xx = sin(pi*x);
else
    xx = pi * cos(pi*x);
end
end

function f = rhs(x) %rhs right hand side
f = pi^2* sin(pi*x(:,1));
end

