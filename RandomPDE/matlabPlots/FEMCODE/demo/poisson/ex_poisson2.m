function ex_poisson2
Nmax = 4;
orders = 3;

[h1,l2,h1s,dof] = deal(zeros(Nmax,orders));

for j=1:orders
	mesh = structured_mesh([0,1,0,1],2,struct('centre',true));
	clear vp;
	vp.elemtype = ['p' num2str(j)];
	vp.a = 1;
	vp.f = @rhs;
	vp.bc = 'zero_dirichlet';
	for i=1:Nmax
		mesh = uniformrefine(mesh);
		[sol, femspace, mesh, vp] = poisson(mesh, vp);
		dof(i,j) = femspace.ndof;
		[h1(i,j), l2(i,j), h1s(i,j)] = ...
			get_H1error_exact(mesh, femspace, sol, @exactu);
	end
end

figure;
loglog(dof(:,1),h1(:,1),'-x', ...
    dof(:,2),h1(:,2),'-x', ...
    dof(:,3),h1(:,3),'-x');
legend('p1','p2','p3');
end


function xx = exactu(x,der)
%exactu exact solution
if nargin<2 || (sum(der)==0)
    xx =sin(pi*x(:,1)) .* sin(pi*x(:,2));
elseif sum(der) == 1 && der(1) == 1
    xx = pi * cos(pi*x(:,1)) .* sin(pi*x(:,2));
else
    xx = pi * sin(pi*x(:,1)) .* cos(pi*x(:,2));
end
end

function f = rhs(x)
%rhs right hand side
f = 2*pi^2* sin(pi*x(:,1)) .* sin(pi*x(:,2));
end

