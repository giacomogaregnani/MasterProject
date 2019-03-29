function [dof, h1, l2, h1s] = ex2_gen

neufun = @(x)(x(:,1) > 1-10^(-5));
dirfun = @(x)(x(:,1) <= 1-10^(-5));
Nmax = 4;
[h1,l2,h1s,dof] = deal(zeros(Nmax,3));
for j=1:3
	mesh = structured_mesh([0,1,0,1],2,struct('centre',true));
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
[H1norm, L2norm, H1seminorm] = get_H1norm(mesh, femspace, sol);
h1 = h1 / H1norm;
l2 = l2 / L2norm;
h1s = h1s / H1seminorm;
figure;
	loglog(dof(:,1),h1(:,1));
hold on;
for i=2:3
	loglog(dof(:,i),h1(:,i));
end
end

function u = exactu(x, der)
	if nargin<2 || (sum(der)==0)
		u = exp(- x(:,1).^2 - x(:,2).^2).*sin(x(:,1).^2 - x(:,2));
	elseif sum(der) == 1 && der(1) == 1
		u = 2.*x(:,1).*exp(- x(:,1).^2 - x(:,2).^2).*cos(x(:,1).^2 - x(:,2)) - 2.*x(:,1).*exp(- x(:,1).^2 - x(:,2).^2).*sin(x(:,1).^2 - x(:,2));
	else
		u = - exp(- x(:,1).^2 - x(:,2).^2).*cos(x(:,1).^2 - x(:,2)) - 2.*x(:,2).*exp(- x(:,1).^2 - x(:,2).^2).*sin(x(:,1).^2 - x(:,2));
	end
end

function f = rhs(x)
	f = zeros(size(x,1),3);
	f(:,1) = -exp(- x(:,1).^2 - x(:,2).^2).*(9.*cos(x(:,1).^2 - x(:,2)) - 8.*exp(x(:,1).^2 + x(:,2).^2) - 19.*sin(x(:,1).^2 - x(:,2)) - 3.*sin(x(:,1) - x(:,2)).*sin(x(:,1).^2 - x(:,2)) + 12.*x(:,2).*cos(x(:,1).^2 - x(:,2)) - x(:,2).*exp(x(:,1).^2 + x(:,2).^2) - 2.*x(:,2).*sin(x(:,1).^2 - x(:,2)) + 2.*cos(x(:,1).^2 - x(:,2)).*cos(x(:,1).*x(:,2)) - exp(x(:,1).^2 + x(:,2).^2).*cos(x(:,1).*x(:,2)) - 2.*sin(x(:,1).^2 - x(:,2)).*cos(x(:,1).*x(:,2)) - 36.*x(:,1).^2.*cos(x(:,1).^2 - x(:,2)) - 4.*x(:,1).^2.*exp(x(:,1).^2 + x(:,2).^2) + 12.*x(:,2).^2.*sin(x(:,1).^2 - x(:,2)) + cos(x(:,1) - x(:,2)).*cos(x(:,1).^2 - x(:,2)) - sin(x(:,1) - x(:,2)).*exp(x(:,1).^2 + x(:,2).^2) - 2.*x(:,1).^2.*x(:,2).*cos(x(:,1).^2 - x(:,2)) + 6.*x(:,1).^2.*x(:,2).*sin(x(:,1).^2 - x(:,2)) - 8.*x(:,1).^2.*cos(x(:,1).^2 - x(:,2)).*cos(x(:,1).*x(:,2)) + x(:,2).*cos(x(:,1) - x(:,2)).*exp(x(:,1).^2 + x(:,2).^2) + 2.*x(:,2).*cos(x(:,1) - x(:,2)).*sin(x(:,1).^2 - x(:,2)) + 4.*x(:,2).*sin(x(:,1) - x(:,2)).*cos(x(:,1).^2 - x(:,2)) - x(:,2).^2.*exp(x(:,1).^2 + x(:,2).^2).*sin(x(:,1).*x(:,2)) - 4.*x(:,1).^2.*x(:,2).^2.*cos(x(:,1).^2 - x(:,2)) + 4.*x(:,1).^2.*x(:,2).^2.*sin(x(:,1).^2 - x(:,2)) + x(:,1).^2.*cos(x(:,1) - x(:,2)).*exp(x(:,1).^2 + x(:,2).^2) + 2.*x(:,1).*x(:,2).*exp(x(:,1).^2 + x(:,2).^2) + 4.*x(:,2).^2.*sin(x(:,1) - x(:,2)).*sin(x(:,1).^2 - x(:,2)) - 2.*x(:,1).*x(:,2).*cos(x(:,1).^2 - x(:,2)).*sin(x(:,1).*x(:,2)) + x(:,1).*x(:,2).*exp(x(:,1).^2 + x(:,2).^2).*sin(x(:,1).*x(:,2)) + 2.*x(:,1).*x(:,2).*sin(x(:,1).^2 - x(:,2)).*sin(x(:,1).*x(:,2)));
	f(:,2) = x(:,1) - x(:,2);
	f(:,3) = x(:,2) + x(:,1).^2;
end

function f = tensor(x,k,l)
	if k==1 && l==1
		f = cos(x(:,1).*x(:,2)) + 5;
	end
	if k==1 && l==2
		f = x(:,1);
	end
	if k==2 && l==1
		f = x(:,1).*x(:,2);
	end
	if k==2 && l==2
		f = sin(x(:,1) - x(:,2)) + 3;
	end
end

function neum = neumann(x, normals)
	if nargin == 1
		error('normals are needed in Neumann B.C. ');
	end
	flux = zeros(size(x,1),2);
	flux(:,1) = - x(:,1).*(x(:,2) + exp(- x(:,1).^2 - x(:,2).^2).*cos(x(:,1).^2 - x(:,2)) + x(:,1).^2 + 2.*x(:,2).*exp(- x(:,1).^2 - x(:,2).^2).*sin(x(:,1).^2 - x(:,2))) - (cos(x(:,1).*x(:,2)) + 5).*(x(:,1) - x(:,2) - 2.*x(:,1).*exp(- x(:,1).^2 - x(:,2).^2).*cos(x(:,1).^2 - x(:,2)) + 2.*x(:,1).*exp(- x(:,1).^2 - x(:,2).^2).*sin(x(:,1).^2 - x(:,2)));
	flux(:,2) = - (sin(x(:,1) - x(:,2)) + 3).*(x(:,2) + exp(- x(:,1).^2 - x(:,2).^2).*cos(x(:,1).^2 - x(:,2)) + x(:,1).^2 + 2.*x(:,2).*exp(- x(:,1).^2 - x(:,2).^2).*sin(x(:,1).^2 - x(:,2))) - x(:,1).*x(:,2).*(x(:,1) - x(:,2) - 2.*x(:,1).*exp(- x(:,1).^2 - x(:,2).^2).*cos(x(:,1).^2 - x(:,2)) + 2.*x(:,1).*exp(- x(:,1).^2 - x(:,2).^2).*sin(x(:,1).^2 - x(:,2)));
	neum = sum(flux .* normals, 2);
end

