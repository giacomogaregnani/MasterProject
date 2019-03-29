function ex_poisson1gen

neufun = @(x)(x(:,1) > 1-10^(-5));
dirfun = @(x)(x(:,1) <= 1-10^(-5));
Nmax = 6;
[h1,l2,h1s,dof] = deal(zeros(Nmax,3));
for j=1:3
	mesh = structured_mesh([0,1],2,struct('centre',true));
	mesh.bdflag{1} = get_bdflag(mesh, dirfun);
	mesh.bdflag{2} = get_bdflag(mesh, neufun);
	clear vp;
	vp.elemtype = ['p' num2str(j)];
	vp.a = @tensor;
	vp.f = @rhs;
	vp.bc{1} = @exactu;
	vp.bc{2} = @neumann;
	for i=1:Nmax
		mesh = uniformrefine(mesh);
		[sol, femspace, mesh, vp] = poisson(mesh, vp);
		dof(i,j) = femspace.ndof;
		[h1(i,j), l2(i,j), h1s(i,j)] = get_H1error_exact(mesh, femspace, sol, @exactu);
	end
end
figure;
loglog(dof(:,1),h1(:,1),dof(:,2),h1(:,2),dof(:,3),h1(:,3));
legend('p1','p2','p3');
end

function u = exactu(x, der)
	if nargin<2 || (sum(der)==0)
		u = exp(x(:,1).^2).*sin(x(:,1));
	elseif sum(der) == 1 && der(1) == 1
		u = exp(x(:,1).^2).*cos(x(:,1)) + 2.*x(:,1).*exp(x(:,1).^2).*sin(x(:,1));
	end
end

function f = rhs(x)
	f = zeros(size(x,1),2);
	f(:,1) = sin(x(:,1)).*(exp(x(:,1).^2).*cos(x(:,1)) - x(:,1) + 2.*x(:,1).*exp(x(:,1).^2).*sin(x(:,1))) - (cos(x(:,1)) + 2).*(exp(x(:,1).^2).*sin(x(:,1)) + 4.*x(:,1).^2.*exp(x(:,1).^2).*sin(x(:,1)) + 4.*x(:,1).*exp(x(:,1).^2).*cos(x(:,1)) - 1);
	f(:,2) = x(:,1);
end

function f = tensor(x,k,l)
		f = cos(x(:,1)) + 2;
end

function neum = neumann(x, normals)
	if nargin == 1
		error('normals are needed in Neumann B.C. ');
	end
	flux = zeros(size(x,1),1);
	flux(:,1) = (cos(x(:,1)) + 2).*(exp(x(:,1).^2).*cos(x(:,1)) - x(:,1) + 2.*x(:,1).*exp(x(:,1).^2).*sin(x(:,1)));
	neum = sum(flux .* normals, 2);
end

