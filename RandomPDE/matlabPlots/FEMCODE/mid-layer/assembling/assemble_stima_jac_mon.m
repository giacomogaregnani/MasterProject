function stima = assemble_stima_jac_mon( mesh, femspace, vp, sol)
% assembles Jacobian matrix for Newton iterations for nonlinear monotone
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
if ~isnumeric(vp.dxiA)
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
                get_dof(mesh, 'all', femspace, j), vp.dxiA(:,i,j), ND, ND);
            continue;
        end

    gradsol=zeros(NT,dim);
    Aij = zeros(NT,1);
    
    for m=1:NQ
        if ~(vp.fully_discrete && strcmp(vp.discrete_type,'tensor'))
            for k=1:dim
                gradsol(:,k) = evalf(mesh, 'all', femspace, lambda(m,:), sol, k);
            end
        end
        for k=1:dim
            for l=1:dim
                if vp.fully_discrete && strcmp(vp.discrete_type,'tensor')
                    vpdxiA = vp.dxiA(:, m, k, l);
                else
                    vpdxiA = vp.dxiA(xloc{m}, gradsol, k, l);
            end
            Aij = Aij + (weight(m) * vpdxiA) .* ...
                evalb(mesh, 'all', lambda(m,:), i, k, ...
                femspace.elemtype) .* ...
                evalb(mesh, 'all', lambda(m,:), j, l, ...
                femspace.elemtype);
        end
        end
    end
    
    Aij = Aij .* mesh.volume;
	stima = stima + sparse(get_dof(mesh, 'all', femspace, i), ...
		get_dof(mesh, 'all', femspace, j), Aij, ND, ND);
    
end
end
end

