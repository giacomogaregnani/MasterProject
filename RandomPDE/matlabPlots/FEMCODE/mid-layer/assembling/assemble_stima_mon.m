function op = assemble_stima_mon( mesh, femspace, vp, sol)
% assembles operator vector for nonlinear monotone elliptic problems

%% INITIALIZATION
dim = size(mesh.elem,2) - 1;
NT = size(mesh.elem,1);
ND = femspace.ndof;
LD = femspace.ldof;
lambda = vp.alambda;
weight = vp.aweight;
NQ = numel(weight);

%% INIT
if ~isnumeric(vp.A)
	xloc = cell([1, NQ]);
	for m=1:numel(weight)
		xloc{m} = get_rc(mesh,[],'all',[],lambda(m,:));
	end
end

%% MATRIX GENERATION
op=zeros(ND, 1);

for i=1:LD
	
    if vp.fully_discrete
		op = op + accumarray(get_dof(mesh, 'all', femspace, i), vp.A(:,i),[ND 1]);
		continue;
    end
    
    gradsol=zeros(NT,dim);
    Ai = zeros(NT,1);
	
    for m=1:NQ
        
        if ~(vp.fully_discrete && strcmp(vp.discrete_type,'tensor'))
        for k=1:dim
            gradsol(:,k) = evalf(mesh, 'all', femspace, lambda(m,:), sol, k);
        end
        end

        for k=1:dim
            if vp.fully_discrete && strcmp(vp.discrete_type,'tensor')
                vpA = vp.A(:, m, k);
            else
                vpA = vp.A(xloc{m}, gradsol, k);
            end
            Ai = Ai + (weight(m) * vpA) .* ...
                evalb(mesh, 'all', lambda(m,:), i, k, ...
                femspace.elemtype);
        end
    end
    
	Ai = Ai .* mesh.volume;
    op = op + accumarray(get_dof(mesh, 'all', femspace, i), Ai,[ND 1]);
end

end

