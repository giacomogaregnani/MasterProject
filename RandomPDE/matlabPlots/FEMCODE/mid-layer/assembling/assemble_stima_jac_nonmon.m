function stima = assemble_stima_jac_nonmon( mesh, femspace, vp, sol)
% assembles Jacobian matrix for Newton iterations for nonlinear nonmonotone
% elliptic problems

%% INITIALIZATION
dim = size(mesh.elem,2) - 1;
NT = size(mesh.elem,1);
ND = femspace.ndof;
LD = femspace.ldof;
lambda = vp.alambda;
weight = vp.aweight;
NQ = numel(weight);

%% INIT
if ~isnumeric(vp.a)
	xloc = cell([1, NQ]);
	for m=1:numel(weight)
		xloc{m} = get_rc(mesh,[],'all',[],lambda(m,:));
	end
end

%% MATRIX GENERATION
stima=sparse(ND, ND);
for i=1:LD
for j=1:LD
	if vp.fully_discrete
		stima = stima + sparse(get_dof(mesh, 'all', femspace, i), ...
			get_dof(mesh, 'all', femspace, j), vp.a(:,i,j), ND, ND);
		continue;
	end
	if vp.symmetric && (j<i)
		continue;
	end
    Aij = zeros(NT,1);
	for m=1:NQ
	for k=1:dim
	for l=1:dim
		if vp.fully_discrete && strcmp(vp.discrete_type,'tensor')
			vpdsa = vp.dsa(:, m, k, l);
		else
			vpdsa = vp.dsa(xloc{m}, evalf(mesh, 'all', femspace, lambda(m,:), sol, 0), k, l);
		end
		Aij = Aij + (weight(m) * vpdsa) .* ...
			evalb(mesh, 'all', lambda(m,:), i, k, ...
			femspace.elemtype) .* ...
			evalf(mesh, 'all', femspace, lambda(m,:), sol, l) .* ...
            evalb(mesh, 'all', lambda(m,:), j, 0, ...
			femspace.elemtype);
	end
	end
	end
	Aij = Aij .* mesh.volume;
	if vp.symmetric && (i<j)
		stima = stima + sparse(get_dof(mesh, 'all', femspace, i), ...
			get_dof(mesh, 'all', femspace, j), Aij, ND, ND);
		stima = stima + sparse(get_dof(mesh, 'all', femspace, j), ...
			get_dof(mesh, 'all', femspace, i), Aij, ND, ND);
	else
		stima = stima + sparse(get_dof(mesh, 'all', femspace, i), ...
			get_dof(mesh, 'all', femspace, j), Aij, ND, ND);
	end
end
end
end

