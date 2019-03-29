function stima = assemble_stima( mesh, femspace, vp )
%ASSEMBLE_STIMA assembles STIFFNESS MATRIX
% A = assemble_stima( mesh, femspace, vp )
%   assembles stiffness matrix A on a mesh, using FE space defined in
%   femspace. The tensor from the divergence form is in the vp structure.

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
        if vp.symmetric && (j<i)
            continue;
        end
        Aij = zeros(NT,1);
        for m=1:NQ
            for k=1:dim
                for l=1:dim
                    if isnumeric(vp.a)
                        if vp.fully_discrete
                            vpa = vp.a(:, k, l, m);
                        elseif numel(vp.a) == 1
                            if k==l,
                                vpa = vp.a;
                            else
                                continue;
                            end
                        else
                            vpa = vp.a(k, l);
                            if vpa == 0
                                continue;
                            end
                        end
                    else
                        vpa = vp.a(xloc{m}, k, l);
                    end
                    Aij = Aij + (weight(m) * vpa) .* ...
                        evalb(mesh, 'all', lambda(m,:), i, k, ...
                        femspace.elemtype) .* ...
                        evalb(mesh, 'all', lambda(m,:), j, l, ...
                        femspace.elemtype);
                end
            end
        end
        Aij = Aij .* mesh.volume;
        stima = stima + sparse(...
            double(get_dof(mesh, 'all', femspace, i)), ...
            double(get_dof(mesh, 'all', femspace, j)), Aij, ND, ND);
        if vp.symmetric && (i<j)
            stima = stima + sparse(...
                double(get_dof(mesh, 'all', femspace, j)), ...
                double(get_dof(mesh, 'all', femspace, i)), Aij, ND, ND);
        end
    end
end
end

