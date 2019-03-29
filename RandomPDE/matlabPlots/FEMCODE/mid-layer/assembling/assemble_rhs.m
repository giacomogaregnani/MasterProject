function rhs = assemble_rhs(mesh, femspace, vp)
%ASSEMBLE_RHS Assembles RIGHT HAND SIDE

if ~isfield(mesh,'volume'), mesh.volume = simplex_volume(mesh); end

%% CONSTANTS
dim = size(mesh.elem,2) - 1;
ND = femspace.ndof;
LD = femspace.ldof;
lambda =  vp.flambda;
weight = vp.fweight;
NQ = numel(weight);

%% ASSEMBLING
rhs=zeros([ND, vp.fdim, vp.NF]);
for j=1:LD
	arhs = 0;
	for i=1:NQ
		%% L2 PART OF FORCE FIELD
        if isnumeric(vp.f)
            if ~vp.fully_discrete_rhs
                eval_fun = vp.f;
            else
                eval_fun = vp.f(:,:,:,i);
            end
        else
            eval_fun = evalf(mesh, 'all', [], lambda(i,:), vp.f);
        end
		arhs = arhs + weight(i) * ...
			bsxfun(@times, eval_fun(:,1:vp.fdim,:), ...
			evalb(mesh,'all',lambda(i,:), j, 0, femspace.elemtype));
		
		%% H-1 PART OF FORCE FIELD (only for Poisson)
		if vp.derivatives
			if vp.fdim > 1
				error('The vector case is not yet implemented');
			end
			for m=1:dim
				for k=1:dim
					if vp.fully_discrete
						vpa = vp.a(:,m,k,i);
					elseif isnumeric(vp.a) % numel(vp.a) == dim^2
						vpa = vp.a(m,k);
					else
						xloc = get_rc(mesh, [],'all',[],lambda(i,:));
						vpa = vp.a(xloc,m,k);
					end
					arhs = arhs + weight(i) * bsxfun(@times, vpa, ...
						bsxfun(@times, eval_fun(:,k+1,:), ...
						evalb(mesh, 'all', lambda(i,:), j, m, ...
						femspace.elemtype)));
				end
			end
		end
	end
	arhs = bsxfun(@times, arhs, mesh.volume);
    dof = get_dof(mesh, 'all', femspace, j);
	for i=1:vp.fdim*vp.NF
		rhs(:,i) = rhs(:,i) + accumarray(dof, arhs(:,i), [ND, 1]);
	end
end
end

