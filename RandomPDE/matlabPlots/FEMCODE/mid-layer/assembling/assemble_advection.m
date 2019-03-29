function stima = assemble_advection( mesh, femspace, vp )
%ASSEMBLE_ADVECTION assembles ADVECTION TERMS
% A = assemble_advection( mesh, femspace, vp )
%   assembles ADVECTION TERMS for system matrix A on a mesh, using FE
%   space defined infemspace. The vector from the advection form is in the
%   vp structure. 

%% INITIALIZATION
dim = size(mesh.elem,2) - 1;
NT = size(mesh.elem,1);
ND = femspace.ndof;
LD = femspace.ldof;
lambda = vp.blambda;
weight = vp.bweight;
NQ = numel(weight);

%% INIT
if ~isnumeric(vp.b)
	xloc = cell([1, NQ]);
	for m=1:numel(weight)
		xloc{m} = get_rc(mesh,[],'all',[],lambda(m,:));
	end
end

%% MATRIX GENERATION
stima=sparse(ND, ND);
for i=1:LD
for j=1:LD
	
    Aij = zeros(NT,1);
	for m=1:NQ
	for k=1:dim
		
        if isnumeric(vp.b)
            if vp.fully_discrete
                vpb = vp.b(:,k,m);
            else
                vpb = vp.b(k);
            end
		else
			vpb = vp.b(xloc{m}, k);
        end
        
		Aij = Aij + (weight(m) * vpb) .* ...
			evalb(mesh, 'all', lambda(m,:), i, k, ...
			femspace.elemtype) .* ...
			evalb(mesh, 'all', lambda(m,:), j, 0, ...
			femspace.elemtype);
	end
	end
	Aij = Aij .* mesh.volume;
	
	stima = stima + sparse(get_dof(mesh, 'all', femspace, j), ...
			get_dof(mesh, 'all', femspace, i), Aij, ND, ND);
end
end
end

