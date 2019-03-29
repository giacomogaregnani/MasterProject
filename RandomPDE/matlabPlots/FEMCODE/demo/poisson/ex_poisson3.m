function ex_poisson3
Nmax = 4;

clear vp sol;
[h1,l2,h1s,dof] = deal(zeros(Nmax,4));

for j=1:3
	mesh = structured_mesh([0,1,0,1,0,1], 1, struct('centre',true));
	
	clear vp;
	vp.elemtype = ['p' num2str(j)];
	vp.a = 1;
	vp.f = @rhs;
	vp.bc = 'zero_dirichlet';
	vp.solver = 'agmg';
	for i=1:Nmax
		[sol, femspace, mesh, vp] = poisson(mesh, vp);
		dof(i,j) = femspace.ndof;
		[h1(i,j), l2(i,j), h1s(i,j)] = ...
			get_H1error_exact(mesh, femspace, sol, @exactu);
		if i<Nmax
			mesh = uniformrefine(mesh);
		end
	end
end

[ H1norm, L2norm, H1seminorm] = get_H1norm(mesh, femspace, sol);
h1 = h1 / H1norm;
l2 = l2 / L2norm;
h2s = h1s / H1seminorm;

figure;
loglog(dof(:,1),h1(:,1),'-x', ...
    dof(:,2),h1(:,2),'-x', ...
    dof(:,3),h1(:,3),'-x');
legend('p1','p2','p3');
end

function xx = exactu(x,der)
%exactu exact solution
if nargin<2 || (sum(der)==0)
    xx =sin(pi*x(:,1)) .* sin(pi*x(:,2)) .* sin(pi*x(:,3));
elseif (sum(der) == 1) && (der(1) == 1)
    xx = pi * cos(pi*x(:,1)) .* sin(pi*x(:,2)) .* sin(pi*x(:,3));
elseif (sum(der) == 1) && (der(2) == 1)
    xx = pi * sin(pi*x(:,1)) .* cos(pi*x(:,2)) .* sin(pi*x(:,3));
elseif (sum(der) == 1) && (der(3) == 1)
    xx = pi * sin(pi*x(:,1)) .* sin(pi*x(:,2)) .* cos(pi*x(:,3));
else
	error('unimplemented derivative');
end
end

function f = rhs(x)
%rhs right hand side
f = 3*pi^2* sin(pi*x(:,1)) .* sin(pi*x(:,2)) .* sin(pi*x(:,3));
end

