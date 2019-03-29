function poisson_example_generator2
syms x y z;
dim = 2;
Nmax = 4;
Omax = 3;
example_name = 'gen';
mesh = 'structured_mesh([0,1,0,1],2,struct(''centre'',true))';

% TENSOR
a = cell(dim,dim);
a{1,1} = 5 + cos(x*y);
a{2,2} = 3 + sin(x-y); 
a{1,2} = x; 
a{2,1} = y*x; 

% SOLUTION
u = exp(-x^2-y^2)*sin(x^2-y); 

% H^-1 force part
g = cell(dim,1);
g{1} = x-y; %x-y;
g{2} = x^2+y; %x^2+y;

%COMPUTE FLUX
flux = cell(dim,1);
flux{1} = simplify(a{1,1}*(diff(u,x) - g{1}) + a{1,2}*(diff(u,y) - g{2})); 
flux{2} = simplify(a{2,1}*(diff(u,x) - g{1}) + a{2,2}*(diff(u,y) - g{2}));

%COMPUTE RHS
f = simplify(-diff(flux{1},x)-diff(flux{2},y)); 

%% OUTPUT
id = fopen(['ex' num2str(dim) '_' example_name '.m'],'w'); % print output
fprintf(id,['function [dof, h1, l2, h1s] = ex' num2str(dim) '_' example_name '\n\n']);
fprintf(id,'neufun = @(x)(x(:,1) > 1-10^(-5));\n');
fprintf(id,'dirfun = @(x)(x(:,1) <= 1-10^(-5));\n');
fprintf(id,['Nmax = ' num2str(Nmax) ';\n']);
fprintf(id,['[h1,l2,h1s,dof] = deal(zeros(Nmax,' num2str(Omax) '));\n']);
fprintf(id,['for j=1:' num2str(Omax) '\n']);
fprintf(id,['\tmesh = ' mesh ';\n']);
fprintf(id,'\tmesh.bdflag{1} = get_bdflag(mesh, dirfun);\n');
fprintf(id,'\tmesh.bdflag{2} = get_bdflag(mesh, neufun);\n');
fprintf(id,'\tclear vp;\n');
fprintf(id,'\tvp.elemtype = [''p'' num2str(j)];\n');
fprintf(id,'\tvp.a = @tensor;\n');
fprintf(id,'\tvp.f = @rhs;\n');
fprintf(id,'\tvp.bc{1} = @exactu;\n');
fprintf(id,'\tvp.bc{2} = @neumann;\n');
fprintf(id,'\tfor i=1:Nmax\n');
fprintf(id,'\t\tmesh = uniformrefine(mesh);\n');
fprintf(id,'\t\t[sol, femspace, mesh, vp] = poisson(mesh, vp);\n');
fprintf(id,'\t\tdof(i,j) = femspace.ndof;\n');
fprintf(id,'\t\t[h1(i,j), l2(i,j), h1s(i,j)] = get_H1error_exact(mesh, femspace, sol, @exactu);\n');
fprintf(id,'\tend\n');
fprintf(id,'end\n');
fprintf(id,'[H1norm, L2norm, H1seminorm] = get_H1norm(mesh, femspace, sol);\n');
fprintf(id,'h1 = h1 / H1norm;\n');
fprintf(id,'l2 = l2 / L2norm;\n');
fprintf(id,'h1s = h1s / H1seminorm;\n');
fprintf(id,'figure;\n');
fprintf(id,'\tloglog(dof(:,1),h1(:,1));\n');
fprintf(id,'hold on;\n');
fprintf(id,['for i=2:' num2str(Omax) '\n']);
fprintf(id,'\tloglog(dof(:,i),h1(:,i));\n');
fprintf(id,'end\n');
fprintf(id,'end\n\n');


%% FUNCTIONS
%SOLUTION (ALSO DIRICHLET B.C.)
fprintf(id,'function u = exactu(x, der)\n');
exactu = substitute(vectorize(char(u)));
dxexactu = substitute(vectorize(char(diff(u,x))));
dyexactu = substitute(vectorize(char(diff(u,y))));
fprintf(id,'\tif nargin<2 || (sum(der)==0)\n');
fprintf(id,'\t\tu = %s;\n',exactu);
fprintf(id,'\telseif sum(der) == 1 && der(1) == 1\n');
fprintf(id,'\t\tu = %s;\n',dxexactu);
fprintf(id,'\telse\n');
fprintf(id,'\t\tu = %s;\n',dyexactu);
fprintf(id,'\tend\n');
fprintf(id,'end\n\n');

%RHS
fprintf(id,'function f = rhs(x)\n');
fprintf(id,'\tf = zeros(size(x,1),%d);\n', dim+1);
rhs = substitute(vectorize(char(f)));
fprintf(id,'\tf(:,1) = %s;\n',rhs);
for i=1:dim
	rhs = substitute(vectorize(char(g{i})));
	fprintf(id,'\tf(:,%d) = %s;\n', i+1, rhs);
end
fprintf(id,'end\n\n');

%TENSOR
fprintf(id,'function f = tensor(x,k,l)\n');
for i=1:dim
	for j=1:dim
		fprintf(id,'\tif k==%d && l==%d\n',i,j);
		aij = substitute(vectorize(char(a{i,j})));
		fprintf(id,'\t\tf = %s;\n', aij);
		fprintf(id,'\tend\n');
	end
end
fprintf(id,'end\n\n');

%FLUX (NEUMANN B.C.)
fprintf(id,'function neum = neumann(x, normals)\n');
fprintf(id,'\tif nargin == 1\n');
fprintf(id,'\t\terror(''normals are needed in Neumann B.C. '');\n');
fprintf(id,'\tend\n');
fprintf(id,'\tflux = zeros(size(x,1),%d);\n', dim);
for j=1:dim
	fl = substitute(vectorize(char(flux{j})));
	fprintf(id,'\tflux(:,%d) = %s;\n', j, fl);
end
fprintf(id,'\tneum = sum(flux .* normals, 2);\n');
fprintf(id,'end\n\n');

%% CLOSE FILE
fclose(id);
end

function ch = vectorize(ch)
ch = strrep(ch,'*','.*');
ch = strrep(ch,'/','./');
ch = strrep(ch,'^','.^');
end

function ch = substitute(ch)
ch = strrep(ch,'exp','EXP');
ch = strrep(ch,'x','x(:,1)');
ch = strrep(ch,'y','x(:,2)');
ch = strrep(ch,'z','x(:,3)');
if isempty(strfind(ch,'x'))
	ch = ['ones(size(x,1),1) * ' ch];
end
ch = strrep(ch,'EXP','exp');
end