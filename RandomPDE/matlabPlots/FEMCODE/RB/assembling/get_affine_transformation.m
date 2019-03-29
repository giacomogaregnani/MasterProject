function [C, G, R] = get_affine_transformation(mesh, ref)
% Takes a mesh that is symbolically defined and computes its affine transf.
%   Detailed explanation goes here

node = mesh.node;
rnode = subsmu(node, ref);
elem  = mesh.elem;

%% CONSTANTS
dim = size(elem,2) - 1;
d = size(rnode,2);
R = size(elem,1);       % number of coarse subdomains

%% AFFINE TRANSFORMATION
% We symbolically precompute the coefficients of affine transformations
% that map the reference domain to the parametrized domains
fprintf('Domain variation - symbolic manipulation\n');

c1 = sym('c1');
c2 = sym('c2');
c3 = sym('c3');
g11 = sym('g11');
g12 = sym('g12');
g13 = sym('g13');
g21 = sym('g21');
g22 = sym('g22');
g23 = sym('g23');
g31 = sym('g31');
g32 = sym('g32');
g33 = sym('g33');

% | x_1 |   | c1 |   | g11 g12 g13 |   | y_1 |
% | x_2 | = | c2 | + | g21 g22 g23 | * | y_2 |
% | x_3 |   | c3 |   | g31 g32 g33 |   | y_3 |
%
% y - point in the reference domain
% x - point in the deformed domain 
% 
% Eqs - euqations
% other notation as in the article

[G,C] = deal(cell(R,1));
for r=1:R
	if (dim == 2) && (d==2)
		Cvec = [c1; c2];
		Gmat = [g11, g12;
			g21, g22];
	elseif (dim == 3) && (d==3)
		Cvec = [c1; c2; c3];
		Gmat = [g11, g12, g13;
			g21, g22, g23;
			g31, g32, g33];
	else
		error 'dimensions undefined';
	end
	
	Eqs = sym(zeros(1,(dim+1)*d));
	for j=1:dim+1
		Eqs(d*(j-1)+1:d*j) = ...
			Cvec + Gmat * transpose(rnode(elem(r,j),:)) ...
			- transpose(node(elem(r,j),:));
	end
	
	if (dim == 2) && (d==2)
		sol = solve(Eqs,c1,c2,g11,g12,g21,g22);
		G{r} = [sol.g11, sol.g12;
			sol.g21, sol.g22];
		C{r} = [sol.c1; sol.c2];
	elseif (dim == 3) && (d==3)
		sol = solve(Eqs,c1,c2,c3,g11,g12,g13,g21,g22,g23,g31,g32,g33);
		G{r} = [...
			sol.g11, sol.g12, sol.g13;
			sol.g21, sol.g22, sol.g23;
			sol.g31, sol.g32, sol.g33];
		C{r} = [sol.c1; sol.c2; sol.c3];
	end
end
R = (1:R)';
end

