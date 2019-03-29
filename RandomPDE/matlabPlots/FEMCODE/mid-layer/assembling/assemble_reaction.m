function stima = assemble_reaction( mesh, femspace, vp )
%ASSEMBLE_REACTION assembles REACTION TERMS
% A = assemble_reaction( mesh, femspace, vp )
%   assembles REACTION TERMS for system matrix A on a mesh, using FE
%   space defined infemspace. The scalar function from the advection form
%   is in the vp structure. 

%% INITIALIZATION
dim = size(mesh.elem,2) - 1;
NT = size(mesh.elem,1);
ND = femspace.ndof;
LD = femspace.ldof;
lambda = vp.clambda;
weight = vp.cweight;
NQ = numel(weight);

%% INIT
if ~isnumeric(vp.c)
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
		if isnumeric(vp.c)
			vpc = vp.c;
			if vpc == 0
				continue;
			end
		else
			vpc = vp.c(xloc{m});
		end
		Aij = Aij + (weight(m) * vpc) .* ...
			evalb(mesh, 'all', lambda(m,:), i, 0, ...
			femspace.elemtype) .* ...
			evalb(mesh, 'all', lambda(m,:), j, 0, ...
			femspace.elemtype);	
	end
	Aij = Aij .* mesh.volume;	
    stima = stima + sparse(...
        get_dof(mesh, 'all', femspace, j), ...
        get_dof(mesh, 'all', femspace, i), Aij, ND, ND);
end
end
end

